import numpy as np
from skimage.filters import gabor_kernel
from scipy.signal import convolve2d
from scipy.ndimage import gaussian_filter
from sklearn.cluster import KMeans

def asm_energy(grid_in, freqs, n_orient):

    ## BUILD THE GABOR FILTER BANK
    orntns = np.deg2rad(np.linspace(0, 157.5, n_orient))  # 8 orntns
    # N = np.min(np.shape(grid_in))
    # max_i = int(np.floor(np.log2(N/8)))
    # i_vals = np.arange(1, max_i + 1)
    # fL = 0.25 - (2 * (i_vals - 0.5)) / N
    # fH = 0.25 + (2 * (i_vals - 0.5)) / N
    # freqs = np.sort(np.array([fL, fH]).flatten())
    # print(freqs)
    # freqs = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5]) # for testing purposes
    # determine the gaussian spread for the gabor filter
    b_set = 0.5 # from (Honarkhah 2011 thesis), and (Pollen & Ronner, 1983)
    lambda_set = 1/freqs
    sig_gauss = lambda_set * (1/np.pi * np.sqrt(np.log(2)/2) * (2**b_set + 1)/(2**b_set - 1))
    # Store kernels in [freq_idx, orientation_idx]
    gabor_kernels = []
    for idx_freq, freq in enumerate(freqs):
        row = []
        for theta in orntns:
            kernel = gabor_kernel(frequency=freq, theta=theta,
                                sigma_x=sig_gauss[idx_freq], sigma_y=sig_gauss[idx_freq])
            row.append(np.real(kernel))
        gabor_kernels.append(row)  # list of rows

    ## APPLY EACH GABOR FILTER IN THE BANK TO THE IMAGE ##
    # define the target image (normalized)
    ti = (grid_in - grid_in.min()) / (grid_in.max() - grid_in.min())
    # ti = grid_in
    # Initialize output: [freq, orientation, H, W]
    filtered_stack = np.zeros((len(freqs), len(orntns), ti.shape[0], ti.shape[1]))
    # Apply all Gabor filters
    for i, freq in enumerate(freqs):
        for j, theta in enumerate(orntns):
            kernel = gabor_kernels[i][j]

            # Step 1: Convolution
            filtered = convolve2d(ti, kernel, mode='same', boundary='symm')

            # Step 2: Nonlinearity (optional, e.g., tanh)
            filtered = np.tanh(0.25 * filtered)

            # Step 3: Gaussian smoothing (optional)
            filtered = gaussian_filter(filtered, sigma=3 * sig_gauss[i])

            # Step 4: Normalize to [0, 1]
            f_min, f_max = filtered.min(), filtered.max()
            if f_max > f_min:
                filtered = (filtered - f_min) / (f_max - f_min)
            else:
                filtered = np.zeros_like(filtered)  # avoid divide-by-zero

            # Store
            filtered_stack[i, j] = filtered

            print(f"Processed Gabor θ={np.rad2deg(theta):.1f}°, f={freq:.3f}")

    # CONVERT THE FILTERED IMAGES TO A LIST FOR FURTHER PROCESSING ##
    filtered_final_stack = np.array(filtered_stack)

    return filtered_final_stack


def asm_cluster(filtered_final_stack, k, spatial_weight):
    ## APPLY K-MEANS CLUSTERING ##
    # Use the correct Gabor stack
    n_freqs, n_orients, H, W = filtered_final_stack.shape
    n_features = n_freqs * n_orients

    # Step 1: Flatten Gabor features to shape (H*W, 64)
    feature_vectors = filtered_final_stack.reshape(n_features, H * W).T  # (H*W, 64)

    # Step 2: Normalize and include spatial coordinates (optional)
    yy, xx = np.meshgrid(np.arange(H), np.arange(W), indexing='ij')  # shape (H, W)

    spatial_weight = spatial_weight
    xx_norm = (xx.flatten() / W) * spatial_weight
    yy_norm = (yy.flatten() / H) * spatial_weight

    # Step 3: Augment features with spatial info
    X_aug = np.concatenate([
        feature_vectors,
        xx_norm[:, np.newaxis],
        yy_norm[:, np.newaxis]
    ], axis=1)  # shape: (H*W, 64+2)

    # Step 4: K-means clustering
    k = k
    seed = np.random.randint(0, 1e6)
    kmeans = KMeans(n_clusters=k, max_iter=300, n_init=1, random_state=seed)
    labels = kmeans.fit_predict(X_aug)

    # Step 5: Reconstruct segmentation image
    segmentation = labels.reshape(H, W)

    return segmentation