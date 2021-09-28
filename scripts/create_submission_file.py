#!/usr/bin/env python
# coding: utf-8
import numpy as np
from hepdata_lib import Submission, Table, Variable, RootFileReader, Uncertainty
import ROOT as rt



submission = Submission()
input_file_dir = '/Users/christinawang/Desktop/Caltech/Research/LLP/SupplementaryMaterials/input_files/'


# Limit plots

limit_table = {}
table_name = ['Figure 3-a', 'Figure 3-b', 'Additional Figure 1']
variable_name = ['c$\\tau$'] + ['95% CL upper limit on B(H $\\rightarrow$ SS)']*6
variable_qualifier = ['c$\\tau$', '-2 $\\sigma$', '-1 $\\sigma$', 'Median expected', '+1 $\\sigma$', '+2 $\\sigma$', 'Observed']
description = 'The 95% CL observed and expected limits on the branching fraction B(H $\\rightarrow$ SS) for 7 GeV mass and $ S \\rightarrow b\\bar{b}$ decay mode.'

for j, decay in enumerate([ 'dddd', '4Tau', 'bbbb']):
    for m in [7, 15, 40, 55]:
        if m == 7 and decay == 'bbbb':continue
        key = decay+str(m)
        limit_table[key] = Table(table_name[j]+' ('+str(m)+' GeV)')
        descrip_temp = description.replace('7 GeV', str(m)+' GeV')
        if decay == 'bbbb':
            limit_table[key].description = descrip_temp
        elif decay == 'dddd':
            limit_table[key].description=descrip_temp.replace('b\\bar{b}', 'd\\bar{d}')
        elif decay == '4Tau':
            limit_table[key].description=descrip_temp.replace('b\\bar{b}', '\\tau^{+} \\tau^{-}')
        limit_table[key].location = table_name[j]
        limit_table[key].keywords["cmenergies"] = ["13000.0"]
        limit_table[key].keywords["observables"] = ["CLS"]
        data = np.loadtxt(input_file_dir+'limits'+decay+'_m'+str(m)+'_hybridNew.txt')
        var = {}
        for i in range(len(variable_name)):#loop over the independent and dependent variables
            if i == 0: 
                var[i] = Variable(variable_name[i], is_independent = True, is_binned = False, units = 'm') #ctau
            else:
                var[i] = Variable(variable_name[i], is_independent = False, is_binned = False, units = '') #limits
                var[i].add_qualifier("Quantile",variable_qualifier[i])


            var[i].values = data[:,i]
            limit_table[key].add_variable(var[i])
        submission.add_table(limit_table[key])



# Region definition & description

region_def = {
'a': 'Region A is defined as 391 cm $< r <$ 695.5 cm and 400 cm $< |z| <$ 671 cm. ',
'b': 'Region B is defined as 671 cm $< |z| <$ 1100 cm, $r <$ 695.5 cm, and $|\eta| <$ 2. ',
}

description = ['The cluster efficiency in bins of hadronic and EM energy in region A. The cluster efficiency is estimated with LLPs decaying to $\\tau^{+} \\tau^{-}$. The sample contains equal fractions of events with LLP mass of 7, 15, 40, and 55 GeV and LLP lifetime of 0.1, 1, 10, and 100m. The first hadronic energy bins correspond to LLPs that decayed leptonically with 0 hadronic energy. The cluster efficiency includes all cluster-level selections described in the paper, except for the jet veto, time cut, and $\Delta\phi$ cut. The full simulation signal yield prediction for samples with various LLP mass between 7 - 55 GeV, lifetime between 0.1 - 100 m, and decay mode to $d\\bar{d}$ and $\\tau^{+} \\tau^{-}$ can be reproduced using this parameterization to within 35% and 20% for region A and B, respectively.',

'The efficiency of $N_{station} > 1$ requirement in bins of hadronic energy in region A.\
 The cluster efficiency is estimated with LLPs decaying to $\\tau^{+} \\tau^{-}$.\
 The sample contains equal fractions of events\
 with LLP mass of 7, 15, 40, and 55 GeV and LLP lifetime of 0.1, 1, 10, and 100m.\
 The first hadronic energy bin corresponds to LLPs that decayed leptonically with 0 hadronic energy.\
 The efficiency is calculated with respect\
 to clusters that pass all cluster-level cuts described in the paper, except for the jet veto, time cut,\
 and $\Delta\phi$ cut.\
 The full simulation signal yield prediction for samples with various LLP mass between 7 - 55 GeV,\
 lifetime between 0.1 - 100 m, and decay mode to $d\\bar{d}$ and $\\tau^{+} \\tau^{-}$ can be reproduced using this\
 parameterization to within 10%.']



#  Cluster Efficiency (Additional Figure 7-8)

files = [ 'cluster_eff_130_correction', 'cutbasedID_eff_130']
names = [ 'Cluster (>130 hits) Efficiency', 'NStation > 1 Efficiency']

for i, f in enumerate(files):
    
    if i == 0:reader = rt.TFile(input_file_dir+f+'.root', 'READ')
    else: reader = RootFileReader(input_file_dir+f+'.root')
    for region in ['a','b']:
        if i == 1 and region == 'a':continue
        if 'cluster' in f:
            hist= reader.Get('h_poly_4Tau_'+region).GetBins()
            histUp= reader.Get('h_poly_4Tau_'+region+'Up').GetBins()
            histDown= reader.Get('h_poly_4Tau_'+region+'Down').GetBins()
            x_edges = []
            y_edges = []
            hist_z = []
            histUp_z = []
            histDown_z = []
            for j in range(len(hist)):
                x_edges.append((hist.At(j).GetXMin(), hist.At(j).GetXMax()))
                y_edges.append((max(hist.At(j).GetYMin(), 0.0), hist.At(j).GetYMax()))
                hist_z.append(hist.At(j).GetContent())
                histUp_z.append(histUp.At(j).GetContent())
                histDown_z.append(histDown.At(j).GetContent())
        else:
            hist = reader.read_hist_1d('h_4Tau_'+region)
            histUp = reader.read_hist_1d('h_4Tau_'+region+'Up')
            histDown = reader.read_hist_1d('h_4Tau_'+region+'Down')
        hadE = Variable('Hadronic Energy', is_independent = True, is_binned = True, units = 'GeV') #hadronic energy            
        EME = Variable('EM Energy', is_independent = True, is_binned = True, units = 'GeV') #EM energy
        eff = Variable('Efficiency', is_independent = False, is_binned = False, units = '') #eff
        eff_unc = Uncertainty('MC statistical uncertainty',is_symmetric=False)
        

        if 'cluster' in f:           
            EME.values = x_edges
            hadE.values = y_edges
            eff.values = hist_z
            eff_unc.values = zip(histUp_z, list(np.array(histDown_z)*-1))
            
        else:
            hadE.values = hist['x_edges']
            eff.values = hist['y']
            eff_unc.values = zip(histUp['y'], list(np.array(histDown['y'])*-1))
        eff.add_uncertainty(eff_unc)    
        if i == 0:table = Table('Additional Figure '+ str(i+7) + region)
        else:table = Table('Additional Figure '+ str(i+7))
        descrip_temp = description[i].replace('in region A. ', 'in region '+region.capitalize()+'. ' + region_def[region])
        if region == 'b': descrip_temp = descrip_temp.replace('35% for region A', '20% for region B')
        
        table.description = descrip_temp
        table.keywords["cmenergies"] = ["13000.0"]
        table.keywords["observables"] = ["EFF"]

        if i == 0: table.location = 'Additional Figure '+ str(i+7) + region
        else: table.location = 'Additional Figure '+ str(i+7)
        table.add_variable(hadE)
        table.add_variable(eff)
        if 'cluster' in f: table.add_variable(EME)


        submission.add_table(table)


# Acceptance plot (Additional Fig 9)

accep_table = Table('Additional Figure 9')
accep_table.description = 'The geometric signal acceptance as a function of $c\\tau$. The acceptance is shown for four different LLP mass hypotheses: 7, 15, 40, and 55 GeV. The acceptance is defined by requiring at least 1 LLP to decay in the region defined as  400 cm $< |z| <$ 1100 cm, $r < $ 695.5 cm, and $|\eta| <$ 2.4.'

accep_table.location = 'Additional Figure 9'
accep_table.keywords["cmenergies"] = ["13000.0"]
accep_table.keywords["observables"] = ["ACC"]

var = {}

for i, m in enumerate([7, 15, 40, 55]):
    data = np.loadtxt(input_file_dir+'acceptance_m'+str(m)+'.txt')

    var[i] = Variable('Acceptance', is_independent = False, is_binned = False, units = '') #limits
    var[i].add_qualifier("LLP mass",str(m)+' GeV')
    var[i].values = data[:,1]
    accep_table.add_variable(var[i])
    
    if m == 55:
        var[0] = Variable('c$\\tau$', is_independent = True, is_binned = False, units = 'm') #ctau
        var[0].values = data[:,0]
        accep_table.add_variable(var[0])

        
    

submission.add_table(accep_table)


# Additional Resources

submission.add_additional_resource('Code for cut-based ID',input_file_dir+'cut_based_id.py', copy_file=True)
submission.add_additional_resource('Signal generator cards',input_file_dir+'generator_cards.tar', copy_file=True)
submission.add_additional_resource('Instructions for Reinterpretation',input_file_dir+'instructions.pdf', copy_file=True)
submission.create_files('submission_file/')

