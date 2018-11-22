# Changelog

## Conventions used for this changelog

 - keep it concise but human readable
 - keep the *UNRELEASED* section up to date with the `develop` branch
 - create a new subsection for each release version
 - each version should have the following information:
   - a release date in the format `YYYY-MM-DD`
   - a list of added new feature
   - a list of changed functionnality of existing features
   - a list of deprecated features (features that will be deleted in a future release)
   - a list of removed feature (previously marked deprecated)
   - a list of bug fixes

----

## *UNRELEASED* (last updated: 2018-11-19)

**Here we collect the list of *added*, *changed*, *deprecated*, *removed* and *fixed* features in preparation for the next release.**

Start of large KKR repository holding *voronoi*, *KKRhost*, *KKRimp*, *KKRsusc*, and *PKKprime* with major refactoring of code structure.


### Added
- masci-tools dependency for kkr_params

### Changed
- kkr calculation retrieves Jij files

### Deprecated
- KKRimporter calculation now also retrieves Jij files

### Removed
- kkr_params deleted from aiida_kkr/tools, instead use masci-tools dependecy

### Fixed
- None

----

