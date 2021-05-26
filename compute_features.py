
feats = ['Z', 'WD', 'batch_Z', 'Entropy', 'Ratio', 'total_PSMs', 'ratio_total_PSMs', 'pep_ratio']

def compute_batch_dist(prey):
    batch_prey_psm_dist = {}
    for batch in prey_batches[prey]:
        psm_ls = batch_prey_psms[batch][prey]
        if all(x == psm_ls[0] for x in psm_ls): psm_ls[0]+=0.5
        psms = psm_ls + [0]*(batch_n_wells[batch] - len(psm_ls))
        batch_prey_psm_dist[(batch, prey)] = (np.mean(psms), np.std(psms))
    return batch_prey_psm_dist

def batch_Z(batch, prey, psm, dist_dict):
    u, se = dist_dict[batch, prey]
    return (psm - u)/se

def bin_feature(feat_vals):
    bins = np.linspace(min(feat_vals)-1e-5, max(feat_vals), 1000)
    return bins