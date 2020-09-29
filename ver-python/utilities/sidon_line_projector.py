import numba as nb
import numpy as np

@nb.njit
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
            
        assert (n_intersect_pts == 2), "Failed to find 2 intersections."


    
    # check if points are not ordered according the direction - flip rows
    if (np.dot(end_pts[1, :] - end_pts[0, :], direction) < 0):
        tmp = end_pts[1].copy()
        end_pts[1] = end_pts[0]
        end_pts[0] = tmp
    
        
    # 2.1 compute intersections x-axes
    x_max = np.max(end_pts[:, 0])
    x_min = np.min(end_pts[:, 0])
    dir_x = end_pts[1, 0] - end_pts[0, 0]
    
    x_axes = np.linspace(-R, R, npixels + 1)[1 : -1]
    x_coords_intersect = x_axes[(x_axes > x_min)*(x_axes < x_max)]
    x_pts_intersect = []
    for x in x_coords_intersect:
        alpha = (x - end_pts[0, 0])  / dir_x
        #assert ((alpha > 0.0) and (alpha < 1.0)), f"Failure to find intersection with line x={x}."
        y = end_pts[0, 1] * (1.0 - alpha) + end_pts[1, 1] * alpha
        x_pts_intersect.append([x, y, alpha]) 
    
    # 2.2 compute intersections y-axes
    y_max = np.max(end_pts[:, 1])
    y_min = np.min(end_pts[:, 1])
    dir_y = (end_pts[1, 1] - end_pts[0, 1])
    
    y_axes = np.linspace(-R, R, npixels + 1)[1 : -1]
    y_coords_intersect = y_axes[(y_axes > y_min)*(y_axes < y_max)]
    y_pts_intersect = []
    for y in y_coords_intersect:
        alpha = (y - end_pts[0, 1])  / dir_y
        #assert ((alpha > 0.0) & (alpha < 1.0))
        x = end_pts[0, 0] * (1.0 - alpha) + end_pts[1, 0] * alpha
        y_pts_intersect.append([x, y, alpha])
    
    # 3. merge ordered arrays of intersection points into one array of ordered intersections points
    #    increasing ordering ordering - column alpha
    
    # order both arrays x_pts_intersect, y_pts_intersect in increasing length with order
    if ((end_pts[1, 0] - end_pts[0, 0]) < 0.0):
        x_pts_intersect.reverse()
        
    if ((end_pts[1, 1] - end_pts[0, 1]) < 0.0):
        y_pts_intersect.reverse()
    
    len_x_intersect = len(x_pts_intersect)
    len_y_intersect = len(y_pts_intersect)
    ind_x = 0
    ind_y = 0
    pts_intersect = []

    # Numba acceleration feature -- append first point in the beginning, pts_insert is not allowed in Numba
    pts_intersect.append([end_pts[0, 0], end_pts[0, 1],  0.0])
    
    
    while (ind_x < len_x_intersect and ind_y < len_y_intersect):
            
        if (x_pts_intersect[ind_x][-1] < y_pts_intersect[ind_y][-1]):
            pts_intersect.append(x_pts_intersect[ind_x])
            ind_x += 1
            continue
            
        if (x_pts_intersect[ind_x][-1] > y_pts_intersect[ind_y][-1]):
            pts_intersect.append(y_pts_intersect[ind_y])
            ind_y += 1
            continue
                
        if (x_pts_intersect[ind_x][-1] == y_pts_intersect[ind_y][-1]):
            pts_intersect.append(x_pts_intersect[ind_x])
            ind_x += 1
            ind_y += 1
            # continue 

    if (ind_x == len_x_intersect):
            
            assert ind_y < len_y_intersect, "Failed to attach the rest of ordered y-array"
            
            for ind in range(ind_y, len_y_intersect):
                pts_intersect.append(y_pts_intersect[ind])
            
    if (ind_y == len_y_intersect):
            assert ind_x < len_x_intersect, "Failed to attach the rest of ordered x-array"
            
            for ind in range(ind_x, len_x_intersect):
                pts_intersect.append(x_pts_intersect[ind])
                
            
    # add last points
    pts_intersect.append([end_pts[1, 0], end_pts[1, 1],  1.0])
    
    # 4. compute the value of the line integral
    ray_int_value = 0.0
    len_intersection = np.linalg.norm(end_pts[0] - end_pts[1])

    for ind in range(len(pts_intersect)-1):
        pix_len = (pts_intersect[ind + 1][2] - pts_intersect[ind][2])
       
        # get indicies of middle points
        mid_x = (pts_intersect[ind + 1][0] + pts_intersect[ind][0]) * 0.5
        mid_y = (pts_intersect[ind + 1][1] + pts_intersect[ind][1]) * 0.5        
        ind_x = int((mid_x + R) / dx) 
        ind_y = npixels - 1 - int((mid_y + R) / dx)

        im_val = 0.0
        if (ind_x >= 0 and ind_x < npixels and ind_y >= 0 and ind_y < npixels):
            im_val = image[ind_y, ind_x]
            
        
        ray_int_value += pix_len * im_val
        
    # computation is over
    return ray_int_value * len_intersection
