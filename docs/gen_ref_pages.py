"""Generate API reference pages from package modules."""

from pathlib import Path

import mkdocs_gen_files


nav = mkdocs_gen_files.Nav()

for path in sorted(Path("repseq").rglob("*.py")):
    module_path = path.with_suffix("")
    doc_path = Path("reference", module_path).with_suffix(".md")
    parts = tuple(module_path.parts)

    if parts[-1] == "__init__":
        doc_path = doc_path.with_name("index.md")
        parts = parts[:-1]
    elif parts[-1] == "__main__":
        continue

    if not parts:
        continue

    nav[parts] = doc_path.relative_to("reference").as_posix()

    with mkdocs_gen_files.open(doc_path, "w") as fd:
        fd.write(f"::: {'.'.join(parts)}\n")
        fd.write("    options:\n")
        fd.write("      show_root_heading: true\n")
        fd.write("      show_source: true\n")

    mkdocs_gen_files.set_edit_path(doc_path, Path("..") / path)

with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
