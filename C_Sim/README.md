
## Building
The binary (`acamar`) is provided for quick testing. However, for testing purposes, the simulator can be built using:
```bash
sh run.sh
```

## Config File
To execute the binary file, a configuration file must be provided for both Acamar and the baseline. The config file should list the clock cycles (extracted from Xilinx HLS synthesis) of different functions used in the simulator. It must contain clock cycles for all functions listed in the `cfg1.txt` file.

## Datasets Format
The `dataset` directory contains datasets in parsed CSR format. The datasets represent 4k x 4k chunks of a larger matrix and must be named as specified in the `dataset` directory.

## Running the Simulator
To run Acamar, execute the `run_sim_host.ipynb` notebook:

```bash
# Launch the notebook
jupyter notebook run_sim_host.ipynb
```
