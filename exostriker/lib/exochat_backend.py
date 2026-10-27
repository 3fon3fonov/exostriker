"""Settings and stateless Responses API requests for the Exo-Striker chat."""
import json
import os
from pathlib import Path

DEFAULTS = dict(model='gpt-6-luna', reasoning_effort='low', verbosity='medium',
                max_output_tokens=8192, api_key='')
MODELS = ('gpt-6-luna', 'gpt-6.1-sol', 'gpt-6-astra')
INSTRUCTIONS = ('You are the Exo-Striker astronomy assistant. Help with exoplanet '
                'analysis, orbital dynamics, statistics, and using Exo-Striker. '
                'Explain units and assumptions when relevant. You only know the '
                'information supplied in this conversation; do not claim to have '
                'inspected or changed the current fit, files, or application. '
                'Give clear, focused answers and distinguish evidence from inference.')


def settings_path():
    folder = Path(__file__).resolve().parent
    developer = folder / 'ES_settings_dev.json'
    from application_paths import user_settings_path
    return developer if developer.is_file() else Path(user_settings_path())


def normalize_settings(raw):
    settings = dict(DEFAULTS)
    settings.update({key: value for key, value in raw.items() if key in DEFAULTS})
    settings['model'] = str(settings['model']).strip() or DEFAULTS['model']
    legacy = settings['model'].startswith('text-davinci-')
    if legacy:
        settings['model'] = DEFAULTS['model']
    if 'max_output_tokens' not in raw and 'max_tokens' in raw and not legacy:
        settings['max_output_tokens'] = int(raw['max_tokens'])
    settings['max_output_tokens'] = int(settings['max_output_tokens'])
    if not 128 <= settings['max_output_tokens'] <= 128000:
        raise ValueError('Output token budget must be between 128 and 128000')
    if settings['reasoning_effort'] not in ('auto', 'low', 'medium', 'high'):
        settings['reasoning_effort'] = DEFAULTS['reasoning_effort']
    if settings['verbosity'] not in ('auto', 'low', 'medium', 'high'):
        settings['verbosity'] = DEFAULTS['verbosity']
    settings['api_key'] = str(settings['api_key']).strip()
    return settings


def read_settings(path=None):
    with Path(path or settings_path()).open() as stream:
        return normalize_settings(json.load(stream).get('openAI', {}))


def save_settings(settings, path=None):
    path = Path(path or settings_path())
    with path.open() as stream:
        document = json.load(stream)
    document['openAI'] = normalize_settings(settings)
    with path.open('w') as stream:
        json.dump(document, stream, indent=2)
        stream.write('\n')


def request_reply(settings, history, message, client_factory=None):
    settings = normalize_settings(settings)
    key = settings['api_key'] or os.environ.get('OPENAI_API_KEY', '').strip()
    if not key:
        raise ValueError('Add an OpenAI API key in Settings or set OPENAI_API_KEY.')
    if client_factory is None:
        from openai import OpenAI
        client_factory = OpenAI
    with client_factory(api_key=key, timeout=120., max_retries=1) as client:
        if not hasattr(client, 'responses'):
            raise RuntimeError('Upgrade the OpenAI SDK: python -m pip install --upgrade "openai>=2.0.0"')
        inputs = list(history) + [dict(role='user', content=message)]
        options = dict(model=settings['model'], input=inputs, instructions=INSTRUCTIONS,
                       max_output_tokens=settings['max_output_tokens'], store=False,
                       include=['reasoning.encrypted_content'])
        if settings['reasoning_effort'] != 'auto':
            options['reasoning'] = dict(effort=settings['reasoning_effort'])
        if settings['verbosity'] != 'auto':
            options['text'] = dict(verbosity=settings['verbosity'])
        response = client.responses.create(**options)
        text = response.output_text
        if not text:
            raise RuntimeError('The model returned no text. If the token budget was exhausted, '
                               'increase it or reduce reasoning effort.')
        output = [item.model_dump(exclude_none=True) for item in response.output]
        notice = ''
        if response.status == 'incomplete':
            notice = 'Response incomplete; increase the output token budget if it was exhausted.'
        return text, inputs + output, notice
