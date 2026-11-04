# Private input manifests

Public `inputs.tsv` files keep the input keys, local paths, assay metadata,
and upstream repository/workflow names. Private GCS addresses are replaced
with `-` in the `gs_uri` column. Public HTTPS references are retained.

Supply private files at the listed `local_subpath` under the input root.
The public manifests cannot download these files from GCS.

Keep a complete private manifest outside the Git checkout. The modular R
scripts accept a separate manifest through these environment variables:

| Figure | Manifest variable | Input root variable |
| --- | --- | --- |
| Figure 1 | `FIG1_MANIFEST` | `FIG1_INPUT` |
| Figure 2 | `FIG2_MANIFEST` | `FIG2_INPUT` |
| Figure 3 | `FIG3_MANIFEST` | `FIG3_INPUT` |
| Figure 4 | `FIG4_MANIFEST` | `FIG4_INPUT` |
| Figure 5 | `FIG5_MANIFEST` | `FIG5_INPUT` |
| Supplement Figure 1 | `SUPPFIG1_MANIFEST` | `SUPPFIG1_INPUT` |

For example:

```sh
export FIG1_MANIFEST=/absolute/path/to/private-manifests/Fig1/inputs.tsv
export FIG1_INPUT=/absolute/path/to/local-inputs
```

The loader uses an existing local file first. A missing file can be downloaded
when the separate manifest supplies an address and the account has access.

The current notebooks read local inputs. Use the input roots in the
[notebook guide](../notebooks/current/README.md). Their recorded execution hashes
describe the original inputs used for the saved run. The manifest hashes from
before URI removal remain in that record; the saved plots were not rerun for
this publication change.

Run `python3 tools/check_public_manifests.py` before publication. GitHub Actions
also checks public files for GCS addresses. The check allows the synthetic
`bucket` host in test scripts and generic scheme descriptions.

Removal from current files does not remove addresses from older Git commits.
