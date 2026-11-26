import sys, glob, os, gemmi, json
root = sys.argv[1]
out = sys.argv[2]
rows = []
for f in sorted(glob.glob(os.path.join(root, 'rcsb_cdr1', '*.cif')) + glob.glob(os.path.join(root, 'pdr5_copy', '*.cif'))):
    doc = gemmi.cif.read(f)
    b = doc.sole_block()
    st = gemmi.read_structure(f)
    st.setup_entities()
    pid = os.path.basename(f)[:-4]
    title = b.find_value('_struct.title')
    res = b.find_value('_refine.ls_d_res_high') or b.find_value('_em_3d_reconstruction.resolution')
    chains = []
    for ch in st[0]:
        poly = [r for r in ch if r.entity_type == gemmi.EntityType.Polymer]
        lig = {}
        for r in ch:
            if r.entity_type != gemmi.EntityType.Polymer and r.name != 'HOH':
                lig[r.name] = lig.get(r.name, 0) + 1
        if poly:
            chains.append(f'{ch.name}:{len(poly)}res({poly[0].seqid.num}-{poly[-1].seqid.num}) lig={lig}')
        elif lig:
            chains.append(f'{ch.name}:nonpoly lig={lig}')
    ref = b.find(['_struct_ref.pdbx_db_accession'])
    accs = [r[0] for r in ref]
    print(pid, '|', title, '|', res, '| UniProt', accs, '| nchains', len(st[0]))
    for c in chains[:8]:
        print('    ', c)
