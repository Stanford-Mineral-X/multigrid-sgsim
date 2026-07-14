from sklearn.cluster import AgglomerativeClustering
from scipy.interpolate import RBFInterpolator
import numpy as np
from .sampling import stratified_sample

def make_trend(fl_xyvc, grid_xyc, smoothing, linespacing):

    # agglomerative clustering
    DistanceThreshold = linespacing/2
    LinkageType = 'average'
    clusteringProgram = AgglomerativeClustering(n_clusters=None, 
                                                #  affinity='euclidean', 
                                                connectivity=None,
                                                linkage=LinkageType,
                                                distance_threshold = DistanceThreshold).fit(fl_xyvc[:,:2])
    labels = clusteringProgram.labels_
    clusterAmount = clusteringProgram.n_clusters_
    stratified_xyv = stratified_sample(fl_xyvc, labels, 1)

    # get x, y, val from stratified sample
    x_stratified = stratified_xyv[:, 0]  # x values
    y_stratified = stratified_xyv[:, 1]  # y values
    val_stratified = stratified_xyv[:, 2]  # values

    # Prepare data for RBF
    data_points = np.column_stack((x_stratified, y_stratified))

    # RBF Interpolation using thin plate spline
    rbf = RBFInterpolator(data_points, val_stratified, kernel='thin_plate_spline', smoothing=smoothing)

    # Interpolate onto the exact same grid points as df_mcgrid
    grdpts = grid_xyc[:, :2]   # (N, 2)
    obspts = fl_xyvc[:, :2]    # (M, 2)
    rbf_values_grdpts = rbf(grdpts)
    rbf_values_obspts = rbf(obspts)

    # add trend values to xyvc arrays on grid points and observation points
    grid_xyct = np.column_stack([grid_xyc, rbf_values_grdpts])
    fl_xyvct = np.column_stack([fl_xyvc, rbf_values_obspts])

    return fl_xyvct, grid_xyct