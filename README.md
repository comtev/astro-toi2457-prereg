# Pre-registered test of a 3.29-day transit signal around TOI-2457

This repository holds a pre-registration, written before the data that will test it were public, together with the code that will run the test. It concerns a possible additional transit signal around TOI-2457 (TIC 117799904), a K dwarf of TESS magnitude 10.87 that already hosts the planet candidate TOI-2457.01 (period 14.2 days). The signal has a period of 3.2911037 days, a depth of about 480 ppm and a duration of about 1.9 hours. If it were a planet transiting this star, it would be roughly 1.3 to 1.4 Earth radii, a figure that is only illustrative because it depends on the adopted stellar radius (0.63 solar radii in the TESS Input Catalog, uncertain by about 10 per cent; catalogue and spectroscopic temperatures for the star also disagree, by up to about 900 K). The signal has not been confirmed, and nothing here shows that it comes from a planet, or from this star rather than a blended neighbour.

The signal was found blind by an amateur multi-sector search of TESS light curves of bright K and M dwarfs, in the QLP light curves of Sectors 32 and 98 (multiple-event statistic 7.6). A separate blind search of the SPOC 2-minute light curves of Sectors 32, 98 and 106 also recovers it (a different pipeline, although the photons are the same). An ephemeris and depth fitted to the QLP data of Sectors 32 and 98 alone predicted the Sector 106 transit times with a 1σ uncertainty of about 3.4 minutes, and a depth of 453 ± 69 ppm; the SPOC light curve of Sector 106 gives 447 ± 103 ppm at those times, roughly 3 to 4σ in one sector depending on processing, which is encouraging but not a confirmation. All three SPOC sectors together give 485 ± 65 ppm, and they fix the cycle count between the sectors (the best alternative count fits worse by Δχ² = 42.5, most of that coming from Sector 106). The signal is still weak, so a further blind test seemed worth setting up before claiming anything.

TESS observes the star again in Sectors 109 (5 to 31 October 2026) and 110 (31 October to 26 November 2026), and the star is on the official 2-minute target lists of both sectors. The pre-registration fixes the predicted times of the 16 transits that fall in those sectors (each with a 1σ uncertainty of about 2.3 minutes), the expected depth, and the rules that decide whether the new data pass the test, fail it, or are inconclusive. Passing would be further evidence that the signal persists; it would not validate a planet or identify the host star. The result will be added to this repository as a new version, whatever it turns out to be.

## Timing of the freeze

The prediction was generated at 02:44 UTC on 11 October 2026 and committed at 02:45:58 UTC (commit 963e1a1 of the analysis code). The pre-registration text was then revised twice the same day, at 09:08:59 UTC (commit c6476de, wording) and at 10:26:06 UTC (commit c7227f3, verdict names, run date and a full test specification), each time changing that one file only; the prediction, the JSON file and the pass and fail criteria did not change, and the text ends with a list of the revisions. Sector 109 observations had started on 5 October 2026, so the first transits had already happened on board, but no Sector 109 or 110 data for this star had been released; `docs/preregistration/mast_inventory_2026-10-11.json` records what MAST returned for the star at 10:28:30 UTC on 11 October 2026.

The commit times above come from a private repository and can only be taken on trust. The independent public timestamps are the first push of this repository to GitHub and the Zenodo record made from its first release (version 1.0.0, [doi:10.5281/zenodo.23297281](https://doi.org/10.5281/zenodo.23297281), created at 11:12:42 UTC), both on 11 October 2026.

## How to cite

Comte, V. 2026, Pre-registered test of a 3.29-day transit signal around TOI-2457 (TIC 117799904) in TESS Sectors 109 and 110, version 1.0.0, Zenodo, [doi:10.5281/zenodo.23297281](https://doi.org/10.5281/zenodo.23297281)

That DOI identifies version 1.0.0, the pre-registration as frozen. The concept DOI [10.5281/zenodo.23297280](https://doi.org/10.5281/zenodo.23297280) stands for all versions of the record and resolves to the latest one, which will include the result once it is added.

## Contents

| Path | Description |
| --- | --- |
| [`docs/preregistration/toi2457_p3.2911.md`](docs/preregistration/toi2457_p3.2911.md) | The pre-registration: the signal, the frozen ephemeris, the predicted transit times, the test rules and a full specification of the test |
| [`docs/preregistration/toi2457_p3.2911.json`](docs/preregistration/toi2457_p3.2911.json) | Machine-readable prediction (ephemeris, alias table, transit template, detrending masks, input files) |
| [`docs/preregistration/toi2457_p3.2911_result_s106_dryrun.json`](docs/preregistration/toi2457_p3.2911_result_s106_dryrun.json) | The test run on Sector 106, which is already part of the prior (so this checks the code, it is not a test of the signal) |
| [`docs/preregistration/mast_inventory_2026-10-11.json`](docs/preregistration/mast_inventory_2026-10-11.json) | What MAST held for the star when the pre-registration was made public |
| [`scripts/prereg_test.py`](scripts/prereg_test.py) | The test |
| [`src/astrohunt/`](src/astrohunt/), [`config/sectors.json`](config/sectors.json) | The parts of the analysis code that the test uses (light-curve loading, detrending, depth fits) |
| [`pyproject.toml`](pyproject.toml), [`uv.lock`](uv.lock) | The environment, with every package pinned to the version that produced the prediction |
| [`scripts/mast_inventory.py`](scripts/mast_inventory.py) | The script that made the MAST inventory (added with this release) |
| [`CODE_MANIFEST.txt`](CODE_MANIFEST.txt) | Commit, tree hash and SHA-256 checksum of every file of the private analysis repository at commit c7227f3 |
| [`SHA256SUMS`](SHA256SUMS) | Checksums of the files in this repository |

The files under `docs/preregistration/` (apart from the MAST inventory), `scripts/prereg_test.py`, `src/astrohunt/` and `config/sectors.json` are byte-for-byte copies of the files of the same names at commit c7227f3, as listed in the manifest. The JSON and dry-run files are identical to those of the first freeze (963e1a1), and so is all the code; only the Markdown text differs. The `git_commit` field inside the JSON (bba21cc) is the parent of 963e1a1, recorded when the prediction was generated. The rest of the analysis code (the search and its other tools) is still private.

## Reproducing the dry run and running the test

With [uv](https://docs.astral.sh/uv/) installed:

```sh
sha256sum -c SHA256SUMS
uv sync --locked
uv run python scripts/prereg_test.py docs/preregistration/toi2457_p3.2911.json 106 --dry-run
```

The last command downloads the SPOC light curves from MAST and rewrites `docs/preregistration/toi2457_p3.2911_result_s106_dryrun.json`; in a clean environment on 11 October 2026 the rewritten file was byte-for-byte identical to the committed one. The test itself will be run, when its data are public, with

```sh
uv run python scripts/prereg_test.py docs/preregistration/toi2457_p3.2911.json 109,110
```

## Caveats

- The depth used as the prior comes from the same sectors in which the signal was found, so it is probably biased upwards to some degree; the consistency test allows for this with a 3σ window.
- A blend with a neighbouring star has not been ruled out. Archival LCO photometry of 1 September 2021, taken while observing TOI-2457.01, covers most of one predicted transit. Under the adopted ephemeris, aperture and error model it disfavours the 24″ neighbour (Gaia G 18.0) and, less strongly (roughly 3 to 4σ, depending on the noise model), TIC 117799899 at 27″ (G 17.3), the neighbour behind TRICERATOPS's nearby eclipsing binary scenario. It does not usefully constrain a source at 32″ (G 19.0), which would need an eclipse of roughly 40 per cent, and a source at 17″ (G 19.2), which would need about 47 per cent, was not measured in it. Fits of the TESS difference images only weakly favour the target over TIC 117799899 (total Δχ² = 1.8) and depend on the assumed point-spread function, so they do not locate the source.
- TRICERATOPS gives a false-positive probability of 0.045 and a nearby false-positive probability of 0.0032, which do not meet the validation thresholds (FPP < 0.015 and NFPP < 10⁻³) of [Giacalone et al. (2021)](https://doi.org/10.3847/1538-3881/abc6af).
- Gaia DR3 reports a RUWE of 1.55 for the star, which might indicate an unresolved companion. No other binarity flag is set in the DR3 fields checked (`ipd_frac_multi_peak` = 0, `non_single_star` = 0), and ExoFOP lists SOAR speckle and Shane adaptive-optics images without a reported companion, but a companion remains possible. Gaia DR4, scheduled for 2 December 2026, may help.

## References

- Brown, T. M., et al. 2013, PASP, 125, 1031, Las Cumbres Observatory Global Telescope Network, [doi:10.1086/673168](https://doi.org/10.1086/673168)
- Gaia Collaboration, Prusti, T., et al. 2016, A&A, 595, A1, The Gaia mission, [doi:10.1051/0004-6361/201629272](https://doi.org/10.1051/0004-6361/201629272)
- Gaia Collaboration, Vallenari, A., et al. 2023, A&A, 674, A1, Gaia Data Release 3: Summary of the content and survey properties, [doi:10.1051/0004-6361/202243940](https://doi.org/10.1051/0004-6361/202243940)
- Giacalone, S., Dressing, C. D., et al. 2021, AJ, 161, 24, Vetting of 384 TESS Objects of Interest with TRICERATOPS and Statistical Validation of 12 Planet Candidates, [doi:10.3847/1538-3881/abc6af](https://doi.org/10.3847/1538-3881/abc6af)
- Hippke, M., David, T. J., et al. 2019, AJ, 158, 143, Wōtan: Comprehensive Time-series Detrending in Python, [doi:10.3847/1538-3881/ab3984](https://doi.org/10.3847/1538-3881/ab3984)
- Huang, C. X., Vanderburg, A., et al. 2020a, RNAAS, 4, 204, Photometry of 10 Million Stars from the First Two Years of TESS Full Frame Images: Part I, [doi:10.3847/2515-5172/abca2e](https://doi.org/10.3847/2515-5172/abca2e)
- Huang, C. X., et al. 2020b, RNAAS, 4, 206, Photometry of 10 Million Stars from the First Two Years of TESS Full Frame Images: Part II, [doi:10.3847/2515-5172/abca2d](https://doi.org/10.3847/2515-5172/abca2d)
- Jenkins, J. M., et al. 2016, Proc. SPIE, 9913, 99133E, The TESS science processing operations center, [doi:10.1117/12.2233418](https://doi.org/10.1117/12.2233418)
- Kreidberg, L. 2015, PASP, 127, 1161, batman: BAsic Transit Model cAlculatioN in Python, [2015PASP..127.1161K](https://ui.adsabs.harvard.edu/abs/2015PASP..127.1161K)
- Kunimoto, M., Huang, C., et al. 2021, RNAAS, 5, 234, Quick-look Pipeline Lightcurves for 9.1 Million Stars Observed over the First Year of the TESS Extended Mission, [doi:10.3847/2515-5172/ac2ef0](https://doi.org/10.3847/2515-5172/ac2ef0)
- Kunimoto, M., Bryson, S., Jaffee, D., et al. 2025, AJ, 170, 280, LEO-Vetter: Fully Automated Flux- and Pixel-level Vetting of TESS Planet Candidates to Support Occurrence Rates, [doi:10.3847/1538-3881/ae070a](https://doi.org/10.3847/1538-3881/ae070a)
- Ricker, G. R., Winn, J. N., Vanderspek, R., et al. 2015, JATIS, 1, 014003, Transiting Exoplanet Survey Satellite, [doi:10.1117/1.JATIS.1.1.014003](https://doi.org/10.1117/1.JATIS.1.1.014003)
- Stassun, K. G., et al. 2019, AJ, 158, 138, The Revised TESS Input Catalog and Candidate Target List, [doi:10.3847/1538-3881/ab3467](https://doi.org/10.3847/1538-3881/ab3467)

## Acknowledgements

This work uses data from the TESS mission, obtained from the Mikulski Archive for Space Telescopes (MAST), including the QLP and SPOC light curves, and data from the European Space Agency mission Gaia, processed by the Gaia Data Processing and Analysis Consortium. It makes use of the [Exoplanet Follow-up Observing Program website](https://exofop.ipac.caltech.edu/tess/) (ExoFOP, [doi:10.26134/ExoFOP5](https://doi.org/10.26134/ExoFOP5)), operated by the California Institute of Technology under contract with NASA. The archival LCO-SSO 1.0 m photometry of 1 September 2021 was obtained and contributed to ExoFOP-TESS by a TFOP observer, to whom we are grateful. An independent review of this pre-registration, with its own re-analysis of the public data, led to the second revision of the text.

## Author

Valentin Comte (independent researcher). The analysis was carried out with the help of Claude, an AI model made by Anthropic; the author checked the results and is responsible for them.

## Licence

The documents and data in this repository are released under the [Creative Commons Attribution 4.0 International licence](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), see [`LICENSE`](LICENSE). The code (`scripts/`, `src/`, `config/`, `pyproject.toml`, `uv.lock`) is released under the [GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html), see [`LICENSE-CODE`](LICENSE-CODE).
