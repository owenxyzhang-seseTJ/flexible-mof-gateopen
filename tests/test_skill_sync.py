from pathlib import Path


def test_standalone_skill_helpers_match_canonical_sources():
    root = Path(__file__).resolve().parents[1]
    skill = root / "skills/flexible-mof-gateopen"
    for copy in (skill / "scripts").glob("*.py"):
        if copy.name in {"quick_validate.py"}:
            continue
        source = root / "src/flexible_mof_gateopen" / copy.name
        if not source.exists():
            source = root / "scripts" / copy.name
        assert copy.read_bytes() == source.read_bytes(), copy.name
    for copy in (skill / "references").glob("*.md"):
        assert copy.read_bytes() == (root / "references" / copy.name).read_bytes()
