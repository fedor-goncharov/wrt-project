def sidon_line_projector(image, phi, shift, radius=1.0):
    """
        Sidon line projector - fast and exact execution algorithm
        for the ray transform along the line

        The code is not yet optimized - works not really fast 
    """
    
    c = np.cos(phi)
    s = np.sin(phi)
    R = radius
    sh = shift
    
    if (sh > R*np.sqrt(2)):
        return 0
    
    # get image size, pixel size
    npixels = image.shape[0]
    dx = 2*R / npixels # pixel's side length 
    
    # set the geometry of the line
    line_center = np.array([[sh * c], [sh * s]])
    direction = np.array([[-s], [c]])
    
    # 1. find intersections of the line with borders [-radius,radius]x[-radius, radius]
    n_intersect_pts = 0
    end_pts = np.zeros((2,2))
    
    # bottom border
    if (np.abs((sh + R*s)) < np.abs(c)*R):
        end_pts[n_intersect_pts, 0] = (sh + R*s)/c
        end_pts[n_intersect_pts, 1] = -R
        n_intersect_pts += 1 
        
    # right border
    if (np.abs((sh - R*c)) < np.abs(s)*R):
        end_pts[n_intersect_pts, 0] = R
        end_pts[n_intersect_pts, 1] = (sh - R*c)/s
        n_intersect_pts += 1 
            
            
    # top border
    if (np.abs((sh - R*s)) < np.abs(c)*R and n_intersect_pts < 2):
        end_pts[n_intersect_pts, 0] = (sh - R*s)/c
        end_pts[n_intersect_pts, 1] = R
        n_intersect_pts += 1 
            
    # left border
    if (np.abs((sh + R*c)) < np.abs(s)*R and n_intersect_pts < 2):
        end_pts[n_intersect_pts, 0] = -R
        end_pts[n_intersect_pts, 1] = (sh + R*c)/s
        n_intersect_pts += 1
            
    assert (n_intersect_pts == 2), f"shift={sh}, phi={phi}, n_intersect={n_intersect_pts}"
    
    
    # check if points are not ordered according the direction - swap rows
    if ( np.dot(end_pts[1, :] - end_pts[0, :], direction) < 0):
        end_pts[[0, 1], :] = end_pts[[1, 0], :] 
    
        
    # 2.1 compute intersections x-axes
    x_max = np.max(end_pts[:, 0])
    x_min = np.min(end_pts[:, 0])
    
    x_axes = np.linspace(-R, R, npixels + 1)[1 : -1]
    x_coords_intersect = x_axes[(x_axes > x_min)*(x_axes < x_max)]
    x_pts_intersect = np.array([], dtype='float64').reshape(0, 3)
    for x in x_coords_intersect:
        alpha = (x - end_pts[0, 0])  / (end_pts[1, 0] - end_pts[0, 0])
        assert ((alpha > 0.0) & (alpha < 1.0))
        y = end_pts[0, 1] * (1.0 - alpha) + end_pts[1, 1] * alpha
        x_pts_intersect = np.append(x_pts_intersect, [[x, y, alpha]], axis=0)
    
    # 2.2 compute intersections y-axes
    y_max = np.max(end_pts[:, 1])
    y_min = np.min(end_pts[:, 1])
    
    y_axes = np.linspace(-R, R, npixels + 1)[1 : -1]
    y_coords_intersect = y_axes[(y_axes > y_min)*(y_axes < y_max)]
    y_pts_intersect = np.array([], dtype='float64').reshape(0, 3)
    for y in y_coords_intersect:
        alpha = (y - end_pts[0, 1])  / (end_pts[1, 1] - end_pts[0, 1])
        assert ((alpha > 0.0) & (alpha < 1.0))
        x = end_pts[0, 0] * (1.0 - alpha) + end_pts[1, 0] * alpha
        y_pts_intersect = np.append(y_pts_intersect, [[x, y, alpha]], axis=0)
    
    # 3. merge ordered arrays into one array of ordered points
    
    # order both arrays x_pts_intersect, y_pts_intersect in increasing length with order
    if ((end_pts[1, 0] - end_pts[0, 0]) < 0.0):
        x_pts_intersect = np.flipud(x_pts_intersect)
        
    if ((end_pts[1, 1] - end_pts[0, 1]) < 0.0):
        y_pts_intersect = np.flipud(y_pts_intersect)
    
    len_x_intersect = x_pts_intersect.shape[0]
    len_y_intersect = y_pts_intersect.shape[0]
    ind_x = 0
    ind_y = 0
    pts_intersect = np.array([], dtype='float64').reshape(0, 3)
    
    if ((len_x_intersect == 0) | (len_y_intersect == 0)):
            pts_intersect = np.concatenate((x_pts_intersect, y_pts_intersect), axis=0)
    else:    
        for _ in range(len_x_intersect + len_y_intersect):
            
            if (ind_x == len_x_intersect):
                pts_intersect = np.append(pts_intersect, y_pts_intersect[ind_y:, :], axis=0)
                break
            
            if (ind_y == len_y_intersect):
                pts_intersect = np.append(pts_intersect, x_pts_intersect[ind_x:, :], axis=0)
                break
            
            if (x_pts_intersect[ind_x, -1] < y_pts_intersect[ind_y, -1]):
                pts_intersect = np.append(pts_intersect, (x_pts_intersect[ind_x, :]).reshape(1,3), axis=0)
                ind_x += 1
                continue
            
            if (x_pts_intersect[ind_x, -1] > y_pts_intersect[ind_y, -1]):
                pts_intersect = np.append(pts_intersect, (y_pts_intersect[ind_y, :]).reshape(1,3), axis=0)
                ind_y += 1
                continue
                
            if (x_pts_intersect[ind_x, -1] == y_pts_intersect[ind_y, -1]):
                pts_intersect = np.append(pts_intersect, (x_pts_intersect[ind_x, :]).reshape(1,3), axis=0)
                ind_x += 1
                ind_y += 1
                continue
                
            
                
    # add first and last points
    start_point = np.array([end_pts[0, 0], end_pts[0, 1],  0]).reshape(1, 3)
    end_point = np.array([end_pts[1, 0], end_pts[1, 1],  1]).reshape(1, 3)
    pts_intersect = np.concatenate((start_point, pts_intersect, end_point), axis=0)
    
    # 4. compute the value of the line integral
    ray_int_value = 0.0
    len_intersection = np.linalg.norm(end_pts[0] - end_pts[1])
    
    for ind in range(pts_intersect.shape[0]-1):
        pix_len = (pts_intersect[ind + 1, 2] - pts_intersect[ind, 2]) * len_intersection
        mid_x = (pts_intersect[ind + 1, 0] + pts_intersect[ind, 0]) * 0.5
        mid_y = (pts_intersect[ind + 1, 1] + pts_intersect[ind, 1]) * 0.5
        ind_x = int((mid_x + R) / dx) 
        ind_y = npixels - int((mid_y + R) / dx) - 1
        
        im_val = 0.0
        if (ind_x >= 0 and ind_x < npixels and ind_y >= 0 and ind_y < npixels):
            im_val = image[ind_x, ind_y]
        
        ray_int_value += pix_len * im_val
        
    # computation is over
    return ray_int_value
