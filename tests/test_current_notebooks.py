"""Check current notebook content, output retention, and portable file paths."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("current", ROOT / "tools/build_current_notebooks.py")
current = importlib.util.module_from_spec(spec)
spec.loader.exec_module(current)


class CurrentNotebooksTest(unittest.TestCase):
    def fixture(self, root):
        folder = root / "notebooks/current"
        folder.mkdir(parents=True)
        (folder / "00_setup.R").write_text('message("setup")\n')
        (folder / "Fig1.R").write_text('# %% Panel D\n#| Rebuild all three assays.\nmessage("plot")\n')
        return {"path": "Fig1/Fig1.ipynb", "script": "notebooks/current/Fig1.R",
                "title": "Current Figure 1", "overview": "Current input values."}

    def test_valid_cells_and_r_kernel(self):
        with tempfile.TemporaryDirectory(prefix="qtl notebook ") as directory:
            root = Path(directory)
            entry = self.fixture(root)
            path, _ = current.build("Fig1", entry, root=root)
            nb = json.loads(path.read_text())
            self.assertEqual(nb["metadata"]["kernelspec"]["name"], "ir")
            self.assertEqual(nb["metadata"]["qtl_reference_date"], "2026-11-02")
            ids = [c["id"] for c in nb["cells"]]
            self.assertEqual(len(ids), len(set(ids)))
            self.assertIn("QTL_PAPER_ROOT", path.read_text())
            self.assertIn("Rebuild all three assays.", path.read_text())

    def test_preserve_outputs_only_when_code_is_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            path, _ = current.build("Fig1", entry, root=root)
            nb = json.loads(path.read_text())
            cell = [c for c in nb["cells"] if c["cell_type"] == "code"][-1]
            cell["outputs"] = [{"output_type": "display_data", "data": {"image/png": "abc"}, "metadata": {}}]
            cell["execution_count"] = 9
            path.write_text(json.dumps(nb))
            current.build("Fig1", entry, root=root)
            saved = json.loads(path.read_text())
            self.assertEqual(saved["cells"][-1]["outputs"], cell["outputs"])
            (root / entry["script"]).write_text('# %% Panel D\nmessage("new plot")\n')
            current.build("Fig1", entry, root=root)
            changed = json.loads(path.read_text())
            self.assertEqual(changed["cells"][-1]["outputs"], [])

    def test_missing_source_does_not_replace_existing_notebook(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            path, _ = current.build("Fig1", entry, root=root)
            before = path.read_bytes()
            (root / entry["script"]).unlink()
            with self.assertRaises(FileNotFoundError):
                current.build("Fig1", entry, root=root)
            self.assertEqual(path.read_bytes(), before)

    def test_changed_setup_clears_downstream_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            path, _ = current.build("Fig1", entry, root=root)
            nb = json.loads(path.read_text())
            for c in nb["cells"]:
                if c["cell_type"] == "code":
                    c["outputs"] = [{"output_type": "stream", "name": "stdout", "text": "old"}]
                    c["execution_count"] = 1
            path.write_text(json.dumps(nb))
            (root / "notebooks/current/00_setup.R").write_text('message("changed setup")\n')
            current.build("Fig1", entry, root=root)
            self.assertTrue(all(not c["outputs"] for c in json.loads(path.read_text())["cells"]
                                if c["cell_type"] == "code"))

    def test_regeneration_keeps_valid_execution_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            path, _ = current.build("Fig1", entry, root=root)
            nb = json.loads(path.read_text())
            provenance = {"R": "4.3.3", "source_inputs_unchanged": True,
                          "input_records": [{"path": "saved.rds", "md5": "abc"}]}
            nb["metadata"]["qtl_execution"] = provenance
            path.write_text(json.dumps(nb))
            current.build("Fig1", entry, root=root)
            self.assertEqual(json.loads(path.read_text())["metadata"].get("qtl_execution"), provenance)
            (root / "notebooks/current/00_setup.R").write_text('message("new setup")\n')
            current.build("Fig1", entry, root=root)
            self.assertNotIn("qtl_execution", json.loads(path.read_text())["metadata"])

    def test_current_index_has_correct_section_numbers(self):
        index = current.read_index(ROOT)
        expected = {f"Fig{i}" for i in range(1, 6)} | {
            "S4", "S11", "S12", "S13", "S14", "S15", "S16", "S20", "S21",
            "S22", "S23", "S24", "S25", "S28", "S30", "S31"}
        self.assertEqual(set(index), expected)
        for entry in index.values():
            self.assertNotIn("/Users/", json.dumps(entry))

    def test_notebook_links_to_actual_preparation_scripts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            entry["preparation"] = [{"key": "fig1_novelty", "script": "preparation/Fig1/prepare_novelty.R",
                                     "description": "Filter source QTL tables and count novel genes."}]
            path, _ = current.build("Fig1", entry, root=root)
            text = path.read_text()
            self.assertIn("prepare_novelty.R", text)
            self.assertIn("QTL_INPUT_ROOT", text)

    def test_includes_show_actual_processing_code(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            helper = root / "notebooks/current/process.R"
            helper.write_text('# %% Filter the source rows\nselected <- dplyr::filter(results, qval < 0.05)\n')
            (root / entry["script"]).write_text('# @include notebooks/current/process.R\n# %% Plot\nmessage("plot")\n')
            path, _ = current.build("Fig1", entry, root=root)
            notebook = json.loads(path.read_text())
            code = "\n".join("".join(c["source"]) for c in notebook["cells"] if c["cell_type"] == "code")
            self.assertIn("dplyr::filter(results, qval < 0.05)", code)
            self.assertNotIn("# @include", code)
            self.assertIn("notebooks/current/process.R", notebook["metadata"]["qtl_source_files"])

    def test_include_change_clears_downstream_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            helper = root / "notebooks/current/process.R"
            helper.write_text('selected <- 1\n')
            (root / entry["script"]).write_text('# @include notebooks/current/process.R\n# %% Plot\nmessage("plot")\n')
            path, _ = current.build("Fig1", entry, root=root)
            nb = json.loads(path.read_text())
            for c in nb["cells"]:
                if c["cell_type"] == "code":
                    c["outputs"] = [{"output_type": "stream", "name": "stdout", "text": "old"}]
                    c["execution_count"] = 1
            nb["metadata"]["qtl_execution"] = {"source_inputs_unchanged": True}
            path.write_text(json.dumps(nb))
            helper.write_text('selected <- 2\n')
            current.build("Fig1", entry, root=root)
            changed = json.loads(path.read_text())
            self.assertNotIn("qtl_execution", changed["metadata"])
            self.assertEqual([c for c in changed["cells"] if c["cell_type"] == "code"][-1]["outputs"], [])

    def test_include_cycles_and_outside_paths_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            script = root / entry["script"]
            script.write_text('# @include notebooks/current/Fig1.R\n')
            with self.assertRaisesRegex(ValueError, "cycle"):
                current.build("Fig1", entry, root=root)
            script.write_text('# @include ../outside.R\n')
            with self.assertRaisesRegex(ValueError, "outside"):
                current.build("Fig1", entry, root=root)

    def test_scope_boilerplate_is_not_shown(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = self.fixture(root)
            entry["limitations"] = ["Internal provenance note for the source guide."]
            path, _ = current.build("Fig1", entry, root=root)
            sources = "".join("".join(c["source"]) for c in json.loads(path.read_text())["cells"]
                              if c["cell_type"] == "markdown")
            self.assertNotIn("## Scope", sources)
            self.assertNotIn("Internal provenance note", sources)


if __name__ == "__main__":
    unittest.main()
