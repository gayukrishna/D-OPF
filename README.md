# D-OPF
This repository contains a toolbox for Distribution Optimal Power Flow (D-OPF) Algorithms developed by Dr. Anamika Dubey's research group (https://eecs.wsu.edu/~adubey/index.html).

## D-OPF Layout

The following is the recommended structure for extracting data from OpenDSS file:

```console
D-OPF
├── D-NET
│   └── network_model.py
│   └── Printed JSON files
├── IEEE-123-Bus
│   └── DSS files
│   └── Printed TXT files
├── LICENSE
└── README.md
```

## Execution

The following procedure will give the JSON files required to model OPF for a distribution network

1. Clone the D-OPF repository
    ```console
    C:\....\> git clone https://github.com/shpoudel/D-OPF
    C:\....\> cd D-OPF
    C:\....\D-OPF>
    ```
1. Run the network_model.py
    ```console
    C:\....\D-OPF> cd D-Net
    C:\....\D-OPF\D-Net>
    C:\....\D-OPF\D-Net> python network_model.py file_type
    
    Note: file_type can be 'JSON' or 'TXT'
    ```
1. Once inside the D-Net folder, you will have the required JSON files and inside the IEEE-123-Bus folder, you will have the required TXT files.
