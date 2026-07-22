# Raw-data layout

Obtain all raw data from the official custodians. Do not commit the files.

Expected layout:

```text
data/raw/
├── faers/
│   ├── archives/
│   │   ├── faers_ascii_2015q1.zip
│   │   └── ...
│   └── manifest/
│       └── download_manifest_sha256.csv
├── canada/
│   └── cvponline_extract_20260331/
└── jader/
    └── extracted/
```

FAERS archives and the SHA-256 manifest can be generated with:

```bash
python scripts/download_faers_ascii.py
```

Canada Vigilance:
https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/canada-vigilance-online-database-data-extract.html

JADER:
https://www.pmda.go.jp/safety/info-services/drugs/adr-info/suspected-adr/0004.html

The analysis used Canada Vigilance through 31 March 2026 and the July 2026 JADER extract. If official filenames change, retain the documented directory names or update the corresponding paths in `src/external_validation.py`.
