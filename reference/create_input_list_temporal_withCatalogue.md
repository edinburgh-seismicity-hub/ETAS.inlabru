# Function to create a default input file for the ETAS Hawkes temporal model where a catalogue is specified in the input file.

Function to create a default input file for the ETAS Hawkes temporal
model where a catalogue is specified in the input file.

## Usage

``` r
create_input_list_temporal_withCatalogue(input_path, num.threads = NULL)
```

## Arguments

- input_path:

  path of the `txt` file containing experiment's information

- num.threads:

  Optional argument for the number of threads to be used by parallel
  processing in inlabru/INLA

## Value

The formatted input.list with the elements required for the temporal
Hawkes model
