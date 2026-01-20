### Catalog Contents

Here, we provide a column-by-column overview of the content of the catalog.
In the following, the available columns are grouped by type.

#### Identifiers

- `ID`: Unique identifier for each object in the catalog.
- `JWST_program_ID`: JWST program ID if the object was observed with JWST (set to nan otherwise).
- `JWST_MSA_ID`: NIRSpec MSA slit ID if the object was observed with JWST (set to nan otherwise).
- `SDSS_plate_mjd_fiberID`: SDSS identifier in the format plate–MJD–fiberID (set to nan for other objects).

#### Sky Coordinates

- `RA`: Right ascension of the object (J2000, in degrees). Set to nan for the stacked spectra.
- `DEC`: Declination of the object (J2000, in degrees). Set to nan for the stacked spectra.

#### Literature References

- `reference`: literature reference, corresponding to the following.
    + Langeroodi+2025 (this work) -> this work
    + Morishita+2024 -> https://ui.adsabs.harvard.edu/abs/2024ApJ...971...43M/abstract
    + Nakajima+2022 -> https://ui.adsabs.harvard.edu/abs/2022ApJS..262....3N/abstract
    + Nakajima+2023 -> https://ui.adsabs.harvard.edu/abs/2023ApJS..269...33N/abstract
    + Revalski+2024 -> https://ui.adsabs.harvard.edu/abs/2024ApJ...966..228R/abstract
    + Sanders+2020 -> https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.1427S/abstract
    + Sanders+2024 -> https://ui.adsabs.harvard.edu/abs/2024ApJ...962...24S/abstract

#### Inferred Properties

- `redshift`: spectroscopic redshift. Set to nan for the stacked spectra.
- `Av`: V-band dust attenuation, in magnitudes. For literature objects where Av measurements are not reported and only attenuation-corrected emission line fluxes are available, this is set to nan.
- `metallicity`: direct-method gas-phase oxygen abundance, 12+log(O/H).
- `metallicity_unc`: 1sigma uncertainty on the direct-method gas-phase oxygen abundance, 12+log(O/H).
- `te(OII)`: Te(OII) electron temperature, in K.
- `te(OII)_unc`: 1sigma uncertainty on the Te(OII) electron temperature, in K.
- `te(OIII)`: Te(OIII) electron temperature, in K.
- `te(OIII)_unc`: 1sigma uncertainty on the Te(OIII) electron temperature, in K.
- `te(OII)_flag`: the method adopted for measuring the Te(OII) electron temperature. "direct" refers to those measured directly using the OII7320,30 flux; "genesis-metallicity" refers to those estimated based on the Te(OIII) value, using the genesis-metallicity non-parametric calibration.

#### Equivalent Widths

- `EW(Hbeta)`: Hbeta rest-frame equivalent width, in Angstroms.
- `EW(Hbeta)_unc`: 1sigma uncertainty on the Hbeta rest-frame equivalent width, in Angstroms.

#### Observed Emission Line Fluxes
(flux units are stored under the `flux_unit` column)

- `observed_O3727,29`: O3727,29 observed flux.
- `observed_O3727,29_unc`: 1sigma uncertainty on the O3727,29 observed flux.
- `observed_Hgamma`: Hgamma observed flux.
- `observed_Hgamma_unc`: 1sigma uncertainty on the Hgamma observed flux.
- `observed_O4363`: O4363 observed flux.
- `observed_O4363_unc`: 1sigma uncertainty on the O4363 observed flux.
- `observed_Hbeta`: Hbeta observed flux.
- `observed_Hbeta_unc`: 1sigma uncertainty on the Hbeta observed flux.
- `observed_O4959`: O4959 observed flux.
- `observed_O4959_unc`: 1sigma uncertainty on the O4959 observed flux.
- `observed_O5007`: O5007 observed flux.
- `observed_O5007_unc`: 1sigma uncertainty on the O5007 observed flux.
- `observed_Halpha`: Halpha observed flux.
- `observed_Halpha_unc`: 1sigma uncertainty on the Halpha observed flux.
- `observed_O7320,30`: O7320,30 observed flux.
- `observed_O7320,30_unc`: 1sigma uncertainty on the O7320,30 observed flux.

- `flux_unit`: units of the reported emission line fluxes.

#### Attenuation-Corrected Line Fluxes
(flux units are stored under the `flux_unit` column)

- `corrected_O3727,29`: O3727,29 attenuation-corrected flux.
- `corrected_O3727,29_unc`: 1sigma uncertainty on the O3727,29 attenuation-corrected flux.
- `corrected_Hgamma`: Hgamma attenuation-corrected flux.
- `corrected_Hgamma_unc`: 1sigma uncertainty on the Hgamma attenuation-corrected flux.
- `corrected_O4363`: O4363 attenuation-corrected flux.
- `corrected_O4363_unc`: 1sigma uncertainty on the O4363 attenuation-corrected flux.
- `corrected_Hbeta`: Hbeta attenuation-corrected flux.
- `corrected_Hbeta_unc`: 1sigma uncertainty on the Hbeta attenuation-corrected flux.
- `corrected_O4959`: O4959 attenuation-corrected flux.
- `corrected_O4959_unc`: 1sigma uncertainty on the O4959 attenuation-corrected flux.
- `corrected_O5007`: O5007 attenuation-corrected flux.
- `corrected_O5007_unc`: 1sigma uncertainty on the O5007 attenuation-corrected flux.
- `corrected_Halpha`: Halpha attenuation-corrected flux.
- `corrected_Halpha_unc`: 1sigma uncertainty on the Halpha attenuation-corrected flux.
- `corrected_O7320,30`: O7320,30 attenuation-corrected flux.
- `corrected_O7320,30_unc`: 1sigma uncertainty on the O7320,30 attenuation-corrected flux.

- `flux_unit`: units of the reported emission line fluxes.
