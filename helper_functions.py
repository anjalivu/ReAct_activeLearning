import numpy as np

def check_inflection(U, conditions, instance):

    has_inflection = False
    for coord_name, coord_value in conditions:
        coord_idx = 0 if coord_name == 'x' else 1
        
        # Find nodes at this boundary and extract deformed coordinates
        deformed_coords = []
        for value in U.values:
            node = instance.nodes[value.nodeLabel - 1]  # nodeLabel is 1-indexed
            if abs(node.coordinates[coord_idx] - coord_value) < 1e-2:
                deformed = [node.coordinates[i] + value.data[i] for i in range(3)]
                deformed_coords.append(deformed)
        
        # Check for inflection if we have enough points
        if len(deformed_coords) > 3:
            deformed_coords = np.array(deformed_coords)
            
            # Sort by z-coordinate
            sort_idx = np.argsort(deformed_coords[:, 2])
            sorted_coords = deformed_coords[sort_idx]
            
            # Get the other coordinate
            z = sorted_coords[:, 2]
            other = sorted_coords[:, 1] if coord_name == 'x' else sorted_coords[:, 0]
            
            # Smooth the coordinates with moving average
            def moving_average(data, window_size):
                return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

            window = 3  # Adjust this (larger = more smoothing)
            other_smooth = moving_average(other, window)
            z_smooth = moving_average(z, window)

            # Calculate derivatives on smoothed data
            d_dz = np.gradient(other_smooth, z_smooth)
            d2_dz2 = np.gradient(d_dz, z_smooth)
            
            # Find inflection points
            sign_changes = np.diff(np.sign(d2_dz2))
            if np.any(sign_changes != 0):
                has_inflection = True
                break

    return has_inflection


def fit_cylinder(displacement_field):
    deformed_coords = []
    for value in displacement_field.values:
        node_label = value.nodeLabel
        
        # Get the instance name and node
        # The value.instance gives us the instance this node belongs to
        instance = value.instance
        
        # Access the node from the instance
        node = instance.nodes[node_label - 1]  # Node labels start at 1, index at 0
        
        # Get original coordinates
        x_orig = node.coordinates[0]
        y_orig = node.coordinates[1]
        z_orig = node.coordinates[2]
        
        # Get displacements
        u1 = value.data[0]
        u2 = value.data[1]
        u3 = value.data[2]
        
        # Calculate deformed coordinates
        x_def = x_orig + u1
        y_def = y_orig + u2
        z_def = z_orig + u3
        
        # Store for cylinder fitting
        deformed_coords.append([x_def, y_def, z_def])
    
    # deformed_coords is assumed to exist and have shape (N, 3)
    points = np.asarray(deformed_coords, dtype=float)
    N = points.shape[0]

    # Axis direction via covariance PCA 
    centroid = points.mean(axis=0)
    X = points - centroid
    cov = np.dot(X.T, X) / N

    eigvals, eigvecs = np.linalg.eigh(cov)
    axis_dir = eigvecs[:, np.argmax(eigvals)]
    axis_dir /= np.linalg.norm(axis_dir)

    # Build orthonormal basis for cross-section plane
    # Choose any vector not parallel to axis
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ref, axis_dir)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])

    e1 = np.cross(axis_dir, ref)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis_dir, e1)
    
    # Project points onto plane normal to axis
    proj = points - np.outer(np.dot(points - centroid, axis_dir), axis_dir)

    x2d = np.dot(proj - centroid, e1)
    y2d = np.dot(proj - centroid, e2)

    # Algebraic circle fit (Taubin-style, NumPy-only)
    x = x2d
    y = y2d

    x_m = x.mean()
    y_m = y.mean()

    u = x - x_m
    v = y - y_m

    Su2 = np.sum(u**2)
    Sv2 = np.sum(v**2)
    Suv = np.sum(u*v)
    Su3 = np.sum(u**3)
    Sv3 = np.sum(v**3)
    Su2v = np.sum(u**2 * v)
    Suv2 = np.sum(u * v**2)

    A = np.array([[Su2, Suv],
                    [Suv, Sv2]])

    B = 0.5 * np.array([Su3 + Suv2,
                        Sv3 + Su2v])

    uc, vc = np.linalg.solve(A, B)

    xc = x_m + uc
    yc = y_m + vc

    radius = np.sqrt((u - uc)**2 + (v - vc)**2).mean()
    
    axis_point = centroid + xc * e1 + yc * e2

    # Goodness-of-fit
    dist = np.sqrt((x - xc)**2 + (y - yc)**2)
    residuals = dist - radius

    rmse = np.sqrt(np.mean(residuals**2))
    mae = np.mean(np.abs(residuals))
    max_error = np.max(np.abs(residuals))

    return radius, dist, residuals