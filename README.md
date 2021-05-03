<!-- PROJECT LOGO -->
<br />
<p align="center">
  <h3 align="center">Genomics Oncology Platform</h3>

  <p align="center">
    project_description
    <br />
    <a href="https://github.com/github_username/repo_name"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name">View Demo</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About the project</a>
      <ul>
        <li><a href="#built-with">Built With </a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

[![Product Name Screen Shot][product-screenshot]](https://example.com)

Genomics Oncology Platform is a graphical user interface software aimed to extract and analyze genomics sequencing data to produce an easily interpretable report. The platform offers tumor purity estimation based on Allele Frequency based Imputation of Tumor Purity (All-FIT; DOI: 10.1093/bioinformatics/btz865), and mutational modeling of variants present in tumor-only sequencing data as well Loss of Heterozygosity (LOH) status based on Loss of Heterozygosity Inference Calculator (LOHGIC; DOI: 10.1200/PO.17.00148).
Paper can be found at: https://www.biorxiv.org/content/10.1101/2021.04.14.439855v1


### Built With

* Python





<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running just clone this repository and run the Platform app in folder dist. 

### Prerequisites

The  user interface allows the user to input a FoundationMedicine xml file or a general input file. All-FIT and LOHGIC are implemented in Python within the user interface.

* Dependencies for the pipeline: the following packages are required for All-FIT and LOHGIC python script
  ```sh
  scipy
  numpy
  textwrap
  os
  xml
  ```
* Dependencies for the GUI:
  ```sh
  tkinter
  pandas
  ```

### Installation

0. Easy run using Anaconda navigator
```sh
Go to: https://www.anaconda.com/products/individual
Download Anaconda navigator based on compatible OS
Launch Spyder
In Spyder's terminal window: cd /path/to/downloaded/git/repository/
Open Platform.py
Click run
```

1. Clone the repo
   ```sh
   git clone https://github.com/njalloul90/Genomics_Oncology_Platform.git
   ```
2. In Mac OS, run the application "Platfrom" in folder dist, OR
3. In windows/linux/mac, the GUI can be run as a python script:
```sh
python Platform.py
```


<!-- USAGE EXAMPLES -->
## Usage

### Using the application

![Alt text](/Screenshots/MainGUI.png?raw=true "Main")


Input file selection:
* General: must be an xlsx file with the required inputs (exact headers: "Sample_ID", "Gene", "VAF", "Depth", "Copy_Number"); the file can include further columns such as pathological and computational purities (exact headers: "Pathological_Purity", "Computational_Purity"), as well as clinical infromation. See /Sample Data/sample_data.xlsx .
* Foundation_xml: must be a FoundationOne CDx xml file with the default tags.


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Your Name - [@nahed_jalloul](https://twitter.com/nahed_jalloul) - email: nj277@cinj.rutgers.edu

Project Link: [https://github.com/njalloul90/Genomics_Oncology_Platform](https://github.com/njalloul90/Genomics_Oncology_Platform)








