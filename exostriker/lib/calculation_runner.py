"""Private subprocess entry point; deliberately avoids importing the GUI."""
import os
import json
from pathlib import Path
import sys
import traceback


def main():
    # Also works from a source checkout; only the installed package's own paths
    # are used when invoked by an installed ProcessWorker.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import job_serialization
    from calculation_job import calculate
    if sys.argv[1] == '--warm':
        folder = Path(sys.argv[2])
        import RV_mod  # Import once, before the user submits a session.
        (folder / 'ready').touch()
        command = sys.stdin.buffer.readline()
        if not command:
            return
        os.chdir(json.loads(command)['cwd'])
        request, response = folder / 'request.dill', folder / 'response.dill'
    else:
        request, response = map(Path, sys.argv[1:])
    try:
        with request.open('rb') as stream:
            fit, operation, options = job_serialization.load(stream)
        result = ('result', calculate(fit, operation, options))
    except BaseException:
        result = ('error', traceback.format_exc())
    temporary = response.with_suffix('.tmp')
    with temporary.open('wb') as stream:
        job_serialization.dump(result, stream)
    os.replace(temporary, response)
    # Keep ownership of any sampling/native descendants until the supervisor has
    # collected the result and cleaned up the entire process tree.
    sys.stdin.buffer.read(1)


if __name__ == '__main__':
    main()
