import re
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).parents[1]
MARKDOWN_FILES = (
    *sorted(ROOT.glob("*.md")),
    *sorted((ROOT / "docs").rglob("*.md")),
    *sorted((ROOT / "literature").rglob("*.md")),
)


def _github_heading_anchors(path: Path) -> set[str]:
    anchors: set[str] = set()
    counts: dict[str, int] = {}
    in_code_fence = False
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("```"):
            in_code_fence = not in_code_fence
            continue
        if in_code_fence or not line.startswith("#"):
            continue
        match = re.match(r"^#{1,6}\s+(.+?)\s*#*$", line)
        if match is None:
            continue
        label = re.sub(r"<[^>]+>", "", match.group(1)).lower()
        base = re.sub(r"[^a-z0-9 _-]", "", label).replace(" ", "-")
        base = re.sub(r"-+", "-", base).strip("-")
        count = counts.get(base, 0)
        counts[base] = count + 1
        anchors.add(base if count == 0 else f"{base}-{count}")
    return anchors


def test_markdown_uses_github_math_delimiters() -> None:
    for path in MARKDOWN_FILES:
        text = path.read_text()
        assert r"\(" not in text, f"{path} uses unsupported inline math delimiters"
        assert r"\)" not in text, f"{path} uses unsupported inline math delimiters"
        assert r"\[" not in text, f"{path} uses unsupported display math delimiters"
        assert r"\]" not in text, f"{path} uses unsupported display math delimiters"
        for token in (r"\,", r"\!", r"\{", r"\}"):
            assert token not in text, (
                f"{path} uses TeX {token}, which GFM consumes before MathJax"
            )
        assert r"\operatorname{" not in text, (
            f"{path} uses operatorname, which GitHub's MathJax policy rejects"
        )
        assert re.search(r"\\\\[^\n]", text) is None, (
            f"{path} has a TeX row separator inside a Markdown line; "
            "put the separator at the end of a source line so GFM preserves it"
        )

        lines = text.splitlines()
        in_code_fence = False
        in_display_math = False
        for index, line in enumerate(lines):
            if line.lstrip().startswith("```"):
                in_code_fence = not in_code_fence
                continue
            if in_code_fence:
                continue
            if line.strip() == "$$":
                assert line == "$$", f"{path}:{index + 1} indents display math"
                if in_display_math:
                    assert index + 1 < len(lines) and not lines[index + 1].strip(), (
                        f"{path}:{index + 1} needs a blank line after display math"
                    )
                else:
                    assert index > 0 and not lines[index - 1].strip(), (
                        f"{path}:{index + 1} needs a blank line before display math"
                    )
                in_display_math = not in_display_math
                continue
            if not in_display_math:
                unescaped = re.sub(r"\\\$", "", line)
                assert unescaped.count("$") % 2 == 0, (
                    f"{path}:{index + 1} has unpaired inline math delimiters"
                )
        assert not in_code_fence, f"{path} has an unclosed code fence"
        assert not in_display_math, f"{path} has an unclosed display-math block"


def test_local_markdown_links_resolve() -> None:
    link_pattern = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
    for source in MARKDOWN_FILES:
        for raw_target in link_pattern.findall(source.read_text()):
            target = raw_target.strip().split(maxsplit=1)[0].strip("<>")
            if target.startswith(("http://", "https://", "mailto:")):
                continue
            path_target, _, anchor = unquote(target).partition("#")
            resolved = (
                (source.parent / path_target).resolve() if path_target else source
            )
            assert resolved.exists(), f"{source} links to missing {path_target}"
            if anchor and resolved.suffix == ".md":
                assert anchor in _github_heading_anchors(resolved), (
                    f"{source} links to missing heading #{anchor} in {resolved}"
                )
