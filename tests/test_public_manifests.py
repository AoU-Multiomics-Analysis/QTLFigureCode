"""Check that published files cannot restore private GCS addresses."""

import importlib.util
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "check_public_manifests",
    Path(__file__).resolve().parents[1] / "tools/check_public_manifests.py",
)
privacy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(privacy)
CLOUD = "gs:" + "//"


class PublicManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def check(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return privacy.check_file(path, self.root)

    def test_local_and_public_https_sources_pass(self):
        self.assertEqual(self.check("Fig1/inputs.tsv", (
            "key\tlocal_subpath\tgs_uri\n"
            "local\tdata/local.tsv\t-\n"
            "public\treference/gencode.gz\thttps://example.org/gencode.gz\n"
        )), [])

    def test_manifest_gcs_values_fail_without_disclosing_address(self):
        address = CLOUD + "private-example-bucket/results/table.tsv"
        issues = self.check("Fig1/inputs.tsv", "key\tgs_uri\nk\t" + address + "\n")
        self.assertEqual(len(issues), 1)
        self.assertIn("Fig1/inputs.tsv:2", issues[0])
        self.assertNotIn(address, issues[0])
        self.assertNotIn("private-example-bucket", issues[0])

    def test_copied_document_address_and_bare_bucket_fail(self):
        for address in (CLOUD + "private-example-bucket/file", CLOUD + "private-example-bucket"):
            with self.subTest(address_type=address.count("/")):
                self.assertTrue(self.check("README.md", "Use `" + address + "`.\n"))

    def test_generic_scheme_and_placeholder_pass(self):
        self.assertEqual(self.check("README.md", "The `gs://` scheme and PRIVATE_GCS_URI.\n"), [])

    def test_only_synthetic_bucket_is_allowed_in_test_scripts(self):
        self.assertEqual(self.check("tests/example.R", CLOUD + "bucket/input.gtf\n"), [])
        self.assertTrue(self.check("tests/example.R", CLOUD + "private-example-bucket/input\n"))
        self.assertTrue(self.check("Fig1/inputs.tsv", CLOUD + "bucket/input\n"))

    def test_binary_files_are_skipped(self):
        path = self.root / "plot.png"
        path.write_bytes(b"\x89PNG\xff\xfe")
        self.assertEqual(privacy.check_file(path, self.root), [])


if __name__ == "__main__":
    unittest.main()
