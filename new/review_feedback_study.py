"""Recheck a completed paired pilot and plot its surface response; no solves."""
import argparse
import json
from pathlib import Path
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    args=parser.parse_args();folder=args.run
    if json.loads((folder/'status.json').read_text())['status']!='completed':
        raise ValueError('Review requires a completed study')
    study=json.loads((folder/'study.json').read_text())
    labels=[x[0] for x in study['cases']]
    if len(labels)!=2:raise ValueError('Expected two cases')
    def load(label,name):return np.loadtxt(folder/label/name,delimiter=',',skiprows=1)
    for name in ('initial_surface.csv','coupling_samples.csv'):
        np.testing.assert_array_equal(load(labels[0],name),load(labels[1],name))
    initial=[load(label,'initial_surface.csv') for label in labels]
    final=[load(label,'final_surface.csv') for label in labels]
    histories=[np.genfromtxt(folder/label/'time_history.csv',delimiter=',',names=True) for label in labels]
    np.testing.assert_array_equal(histories[0]['time_s'],histories[1]['time_s'])
    np.testing.assert_array_equal(final[0][:,0],final[1][:,0])
    delta=final[0][:,1]-final[1][:,1]
    comparison=json.loads((folder/'comparison.json').read_text())
    np.testing.assert_allclose(np.sqrt(np.mean(delta**2)),comparison['surface_rms_difference_nm'],atol=1e-14,rtol=1e-10)
    report={'source':str(folder),'cases':labels,'initial_surface_and_energy_equal':True,
        'surface_rms_difference_nm':float(np.sqrt(np.mean(delta**2))),
        'surface_max_difference_nm':float(np.max(np.abs(delta))),
        'maximum_difference_radius_nm':float(final[0][np.argmax(np.abs(delta)),0]),
        'evolution':{label:{'center_change_nm':float(b[0,1]-a[0,1]),
            'edge_change_nm':float(b[-1,1]-a[-1,1]),
            'max_absolute_phase_mass_drift':float(np.max(np.abs(h['relative_mass_drift'])))}
            for label,a,b,h in zip(labels,initial,final,histories)},
        'scope':'Verification and plotting of saved results; no simulation or physical validation'}
    (folder/'review.json').write_text(json.dumps(report,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(7,9))
    for label,a,b in zip(labels,initial,final):axes[0].plot(b[:,0],b[:,1]-a[:,1],label=label)
    axes[0].set(xlabel='Reference radius (nm)',ylabel='Contour height change (nm)');axes[0].legend()
    axes[1].plot(final[0][:,0],delta)
    axes[1].set(xlabel='Reference radius (nm)',ylabel='Surface difference (nm)',title=f'{labels[0]} minus {labels[1]}')
    h0,h1=histories
    axes[2].plot(h0['time_s'],h0['edge_z_phase_nm']-h1['edge_z_phase_nm'],label='Edge')
    axes[2].plot(h0['time_s'],h0['center_z_phase_nm']-h1['center_z_phase_nm'],label='Center')
    axes[2].set(xlabel='Time (s)',ylabel='Surface difference (nm)');axes[2].legend()
    fig.tight_layout();fig.savefig(folder/'feedback_review.png',dpi=170);plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
