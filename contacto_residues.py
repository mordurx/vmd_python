#%%
def residue_for_atoms_id(atom_list, topology):
    return set([topology.atom(a).residue.index for a in atom_list])
def residue_for_atoms_name(atom_list, topology):
    vect=[]
    for a in atom_list:
        vect.append(topology.residue(a))
    return vect
def residue_for_atoms_original(atom_list, topology):
    vect=[]
    index=90
    for a in atom_list:
        if (index==354):
            print ("corte")
            index=476
        if (index==701):
            print ("corte2")
            index=1148    
        if (index==1432):
            print ("corte3")
            index=1478

        resid=str(topology.residue(a).code+str(index))
        index=index+1
        
        vect.append(resid)
    return vect    
def segment_for_residue(topology):
    vect=[]
   
    R1 = topology.select("segname R1")
    R2 = topology.select("segname R2")
    R3 = topology.select("segname R3")
    R4 = topology.select("segname R4")
    for a in R1:
        vect.append(topology.residue(a))
    return vect

import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
from contact_map import ContactMap, ContactFrequency, ContactDifference
pdb='poses/snx_chanel'
traj = md.load_pdb(pdb+'.pdb')
print(traj)
topology = traj.topology

tox = topology.select("segname TOX")
cav = topology.select("segname R1 R2 R3 R4")

frame_contacts= ContactMap(traj[0],query=tox,haystack=cav,cutoff=0.35)
#print (frame_contacts.residue_contacts.df)
df=frame_contacts.residue_contacts.df

(fig, ax) = frame_contacts.residue_contacts.plot(cmap='seismic', vmin=-1, vmax=1)

tox_residues_id = residue_for_atoms_id(tox, traj.topology)
cav_residues_id = residue_for_atoms_id(cav, traj.topology)
tox_residues = residue_for_atoms_name(tox_residues_id, traj.topology)
cav_residues = residue_for_atoms_name(cav_residues_id, traj.topology)
cav_residues_ori = residue_for_atoms_original(cav_residues_id, traj.topology)
ax.set_xlim(min(cav_residues_id), max(cav_residues_id) + 1)
ax.set_ylim(min(tox_residues_id), max(tox_residues_id) + 1)
#segment_for_residue(topology)
#myDataFrame.set_index('column_name')

df=df.iloc[min(tox_residues_id):max(tox_residues_id) + 1, min(cav_residues_id):max(cav_residues_id) + 1]
df.columns = [cav_residues,cav_residues,cav_residues_ori]
columns = list(tox_residues)
#p = pd.reindex(columns=columns)
df.insert(0,'tox', columns)
df=df.set_index('tox')
df=df.T
#columns = pd.MultiIndex.from_tuples([('speed', 'max'),
#...                                      ('species', 'type')])
#indexedData = data.set_index('Athlete')
#df = pd.DataFrame(l[min(tox_residues):max(tox_residues) + 1][min(cav_residues): max(cav_residues)+ 1],columns = cav_residues)


# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter(pdb+'.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name='Sheet1')
# Get the xlsxwriter workbook and worksheet objects.
workbook  = writer.book
worksheet = writer.sheets['Sheet1']
start_row = 1
start_col = 1
end_row = max(cav_residues_id)+ 1
end_cold = max(cav_residues_id)+ 1

# Add a format. Light red fill with dark red text.
format1 = workbook.add_format({'bg_color': '#FFC7CE',
                               'font_color': '#9C0006'})
# Apply a conditional format to the cell range.
worksheet.conditional_format(start_row, start_col, end_row, end_cold,
                             {'type': 'no_blanks',
                              
                              
                              'format':   format1})

# Close the Pandas Excel writer and output the Excel file.
writer.save()


plt.xlabel("Residue")
plt.ylabel("Residue")
plt.savefig('atom_contacts.pdf')

plt.show()






# %%

