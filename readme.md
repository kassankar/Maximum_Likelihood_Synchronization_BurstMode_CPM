
# Maximum Likelihood Synchronization of Burst-Mode CPM

A Matlab code for the implementation of the synchronization method introduced by __Hosseini & Perrins__  for different types of CPM signals (GMSK-RECT-RC....). The details of the implementation can be obtained from this two papers: _"Maximum Likelihood Synchronization of Burst-Mode CPM"_ and _"Timing, Carrier, and Frame Synchronization of Burst-Mode CPM"_


# How to cite this Work ?
If you use this package for your own work, please consider citing it with [this piece of BibTeX](DATAAIDED_SSB.bib):

```bibtex
@INPROCEEDINGS{9615878,
  author={Kassan, Karim and Farès, Haïfa and Louët, Yves and Glattli, Christian},
  booktitle={2021 International Symposium on Networks, Computers and Communications (ISNCC)}, 
  title={Training Sequence Design For Burst Mode Single Side Band CPM Synchronization}, 
  year={2021},
  volume={},
  number={},
  pages={1-5},
  doi={10.1109/ISNCC52172.2021.9615878}}
```
# How to run the code ?
1. Make sure that you have a compatibel version of python (this code was tested on python 3.8)
2. Download (clone) the files from the repository.
3. Open the file called _MLS_BurstMode_CPM.py_
4. Select the section called __Variables ini__
5. Select the type of pulses by changing the variable `pulse` number
	* `1` is for Lorentzian
	* `2` is for GMSK
	* `3` is for Raised Cosine
	* `4` is for Rectangular
6. Change the pulse length by changing the variable `L`
7. Select the sampling frequency by changing the variable `os` ( default is `2`)
8. Select M-ary by changing the `M` ( e.g M=2 for Binary)
9. Select the modulation index by changing the `h`.
10. Select the section called __Training sequence Phase and Signal__
    * Select the data added preamble length `L0`.
12. Select the section called __Add channel Noise and offset__
    * select the required SNR values from `Snr`.
    * select the required number of iterations `itt` (default is `10**5`).
## License
© 2020-2021 Karim Kassan

