"""Package resources and writable per-user application settings."""
from pathlib import Path
import shutil


def package_resource(relative_path):
    """Return an absolute path within the installed lib directory."""
    return str(Path(__file__).resolve().parent / relative_path)


def user_settings_path():
    """Seed user settings from the packaged defaults on first use."""
    directory = Path.home() / '.exostriker'
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / 'ES_settings.json'
    if not destination.exists():
        shutil.copyfile(package_resource('ES_settings.json'), destination)
    return str(destination)
