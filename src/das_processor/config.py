"""Configuration helpers for das_processor.

The package ships two editable configuration templates, ``defaults.py`` and
``sources.py``, in ``das_processor/templates/``. ``init_config`` copies them
into the user's working folder, where they can be edited and imported with::

    from sources import sources
    from defaults import channel_ranges, path_params, ...
"""

import importlib
import importlib.util
from importlib.resources import files
from pathlib import Path

TEMPLATES = ("defaults.py", "sources.py")


def init_config(dest=".", overwrite=False):
    """Copy the configuration templates into ``dest``.

    Existing files are never overwritten unless ``overwrite=True``.

    Parameters
    ----------
    dest : str or Path
        Target folder (created if needed). Default: current directory.
    overwrite : bool
        Replace files that already exist.

    Returns
    -------
    list[Path]
        Paths of the files that were created.
    """
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    template_dir = files(__package__).joinpath("templates")

    created = []
    for name in TEMPLATES:
        target = dest / name
        if target.exists() and not overwrite:
            print(f"Skipping {target} (already exists)")
            continue
        target.write_bytes(template_dir.joinpath(name).read_bytes())
        print(f"Created {target}")
        created.append(target)

    # Make freshly created files importable in the same Python session
    importlib.invalidate_caches()
    return created


def load_module(path, name):
    """Load a Python file as a module, without adding it to ``sys.modules``."""
    path = Path(path).resolve()
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `das-processor init` or `init_config()` "
            "to create the template files."
        )
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_config(config_dir="."):
    """Load ``defaults.py`` and ``sources.py`` from ``config_dir``.

    Returns
    -------
    tuple
        ``(defaults_module, sources_module)``
    """
    config_dir = Path(config_dir)
    defaults = load_module(config_dir / "defaults.py", "user_defaults")
    sources = load_module(config_dir / "sources.py", "user_sources")
    return defaults, sources
