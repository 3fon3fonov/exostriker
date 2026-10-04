# Exo-Striker AI Chat

Exo-Striker 0.96 uses the OpenAI Responses API instead of the legacy Davinci
completion endpoint. Upgrade the SDK in the Python environment used to run
Exo-Striker. The chatbot is optional; a normal installation does not install
the OpenAI SDK. Enable it with the `chat` extra:

```sh
python -m pip install --upgrade 'exostriker[chat]'
```

From a Git checkout, use `python -m pip install '.[chat]'`. You can also
install just the SDK with `python -m pip install --upgrade 'openai>=2.0.0'`.

Open the **Exo-Striker AI Chat** tab and click **Settings**. Supply your OpenAI
API key, or leave it blank to use `OPENAI_API_KEY`. The key field is masked;
keys entered there are saved in the existing local settings JSON. Using the
environment variable avoids storing a key in that file.

The model picker offers GPT-6 Luna (default, fast), GPT-6.1 Sol (balanced), and
GPT-6 Astra (complex work). Model access depends on your API account. You can
enter another API model ID without editing the GUI code. For a custom model
that does not support reasoning or response-detail settings, select **auto**
for those fields to omit the corresponding request parameters.

Suggested starting settings are low reasoning effort, medium response detail,
and an output token budget of 8192. The budget includes reasoning and visible
answer tokens; it is not a word count. More reasoning can increase latency and
token use. Increase the budget if a response is incomplete. The old temperature,
top-p and repetition-penalty controls are removed; they are not a portable
request contract for current reasoning models.

Requests run in a worker thread, keeping the GUI responsive. Successful replies
carry the conversation forward, including encrypted reasoning output items.
**New chat** clears the local history. Only typed chat messages and prior replies
are sent: current fit parameters, files and plots are not attached automatically.
The application requests `store=False`; this is not a claim of zero retention
under all OpenAI API policies. No code generated in a reply is executed.

Old Davinci settings migrate to the current defaults when read, while retaining
an existing API key. Saving replaces the old OpenAI options and preserves the
other settings sections.

Official references:

- [Current model guidance](https://developers.openai.com/api/docs/guides/latest-model)
- [Migrating to Responses](https://developers.openai.com/api/docs/guides/migrate-to-responses)
- [Conversation state](https://developers.openai.com/api/docs/guides/conversation-state)
