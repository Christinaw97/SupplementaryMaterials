"""

This file provides a parameterized implementation of the cut-based ID described in the 
CMS analysis "Search for long-lived particles decaying in the endcap muon system in proton-proton 
collisions at sqrt{s} = 13 TeV". The provided parameterization can reproduce the cut-based ID efficiency from full simulation
to within 10% for various LLP masses (7-55 GeV), lifetimes (0.1 - 100 m), and decay modes (a pair of d quarks or a pair of tau leptons) considered.
A detailed instruction on the usage of the code is provided in the 
HEPData record: https://doi.org/10.17182/hepdata.104408

"""


import numpy as np
import math
import ROOT as rt
"""

The AvgStation takes arrays of LLP decay position z and r in cm, outputs the avgStation prediction. 
More details in the HEPData entry.

"""
 
def AvgStation(z, r): 
    station = np.copy(z)
    station[np.abs(z)<632] = 1
    station[np.logical_and(np.abs(z)<724, np.abs(r)>275)] = 1
    station[np.logical_and(station>1, np.abs(z)<850)] = 2
    station[np.logical_and(np.abs(z)>=850,np.abs(z)<970)] = 3
    station[np.abs(z)>=970] = 4
    return station

"""

The region function outputs the LLP decay region as a string of 'a', 'b', or 'NA', depending on the LLP decay position z and r in cm and the LLP eta

"""
def region (z_list, r_list, eta_list): 
	region_list = []
	for i in range(len(z_list)):
		z = z_list[i]
		r = r_list[i]
		eta = eta_list[i]
		if r > 391 and r < 695.5 and abs(z) > 400 and abs(z) < 671: region_list.append('a')
		elif abs(z) > 671 and abs(z) < 1100 and r < 695.5  and abs(eta) < 2: region_list.append('b')
		else: region_list.append('NA')
	return np.array(region_list)

"""

input: list of LLP decay position in z, r, LLP hadronic energy, LLP eta, and 
a histogram of parametrization of NStation>1 efficiency with respect to the LLP hadronic energy (Additional Figure 8 in HEPData record)
input type: z, r, hadE, and eta are numpy arrays and eff_hist is TH1F
output: list of probability that the event passes cut-based ID


"""

def cut_based_id(z, r, hadE, eta, eff_hist_b):
	avgStation = AvgStation(z, r)
	regions = region(z, r, eta)
	eta_cut = np.copy(avgStation)
	eta_cut[eta_cut==1]  = 1.8 #implicitly 1.1
	eta_cut[eta_cut==2]  = 1.6
	eta_cut[eta_cut==3]  = 1.6
	eta_cut[eta_cut==4]  = 1.8
	weight=[]
	for j in range(len(hadE)):  
		if regions[j] == 'NA': 
			weight.append(0.0)
		elif regions[j] == 'a':
			weight.append(1.0)
		else: 
			x = eff_hist_b.GetXaxis().FindFixBin(hadE[j])
			weight.append(eff_hist_b.GetBinContent(x))
	weight = np.array(weight)

	#apply eta<1.9 if NStation>1, else apply different eta cut as defined above
	weight = weight*(np.abs(eta)<1.9)+(1-weight)*(np.abs(eta)<eta_cut)  
	weight[regions[j]=='NA'] = 0.0
	return weight

if __name__ == "__main__":
	"""
	dummy IO of input lists and efficiency parameterization (Additional Figure 8).
	Path to ROOT file and input lists z, r, hadE, and eta need to be replaced. 

	"""

	rootFile = '/PATH/Additional_Figure_8.root'
	directory_name = 'Additional Figure 8'
	hist_name = 'Hist1D_y1'
	z = [ 865.03613,  535.4038 ,  473.5255 , 1065.5292 ]
	r = [329.1767 , 660.40015, 296.6402 , 313.69202]
	hadE = [81.33661  , 43.15993  ,  3.8601131, 67.47514 ]
	eta = [1.6940517 , 0.7383448 , 1.2472409 , 1.9405304]
	file = rt.TFile(rootFile, 'READ')
	eff_hist_b = file.Get('Additional Figure 8').Get('Hist1D_y1')
	assert(len(z) == len(r) == len(hadE) == len(eta))


	# main function execution
	weights = cut_based_id(z, r, hadE, eta, eff_hist_b)
	print("cut-based ID efficiency: ", np.sum(weights)/len(weights))

	file.Close()
