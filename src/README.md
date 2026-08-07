These are the source files for the package.
They contain the following content

* structs.jl: Defines the FusionRing type and its constructor
* creation.jl : Functions for constructing fusion rings of various kinds.
* FusionRings.jl : Main file. Loads fusion ring data and initializes symbolic numbers. At the moment it contains only the essentials. Keep as simple as possible.
* operations.jl: Functions to manipulate fusion rings.
* properties.jl: Functions that query and/or calculate properties of fusion rings
* general_functions.jl: Helper functions that don't belong in any of the other categories
* import_data.jl: Functions to import and export FusionRing objects and symbolic numbers
* formatting_and_printing: Functions to format and print data
* categorifiability_criteria: functions to check whether a fusion ring can't be categorified to a fusion category with certain properties
