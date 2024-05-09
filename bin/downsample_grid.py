
import numpy as np

def downsamp_reg_grid(old_mat, old_x_edges, old_y_edges, new_y_edges, new_x_edges, no_data_val=-99999,
                      da=None):
    # Downsample a regular grid by integration to another regular grid.  The 'old' grid has axis grid
    # edges defined by old_x_edges and old_y_edges.  The 'new' (output) grid is defined by new_x_edges
    # and new_y_edges.  The function excludes no_data_val pixels from the integration, and furthermore
    # sets new pixels made up from > 0.5 portions of no_data_val old-pixels to no_data_val.
    # In the event that the pixel-areas are different than dx*dy from the rectilinear axes, the user
    # may input pixel areas of the old grid using the 'da' input.  da should be a matrix that matches
    # the dimensions of old_mat.

    # determine overlap weights for each row and column of new grid
    #   include area-weighting in row associations
    new_y_n = len(new_y_edges) - 1
    old_y_n = len(old_y_edges) - 1
    old_y_widths = np.diff(old_y_edges)
    row_weight_mat = np.zeros((new_y_n, old_y_n), dtype=float)
    row_da_weight  = np.zeros((new_y_n, old_y_n), dtype=float)
    for new_y_index in range(new_y_n):
        # determine linear row-weighting of original pixels to new pixels
        temp_edges = new_y_edges[new_y_index:(new_y_index+2)]
        pixel_portions = pixel_portion_overlap1D(temp_edges, old_y_edges)
        bin_indices = np.where(pixel_portions > 0.)
        bin_weights = pixel_portions[bin_indices]
        row_da_weight[new_y_index, bin_indices] = bin_weights
        # also weight by pixel width(height).
        area_weights = old_y_widths[bin_indices]
        bin_weights = bin_weights*area_weights
        # normalize
        bin_weights = bin_weights/bin_weights.sum()
        # store indices and weights for each row
        row_weight_mat[new_y_index, bin_indices] = bin_weights

    # repeat for columns
    new_x_n = len(new_x_edges) - 1
    old_x_n = len(old_x_edges) - 1
    old_x_widths = np.diff(old_x_edges)
    column_weight_mat = np.zeros((old_x_n, new_x_n), dtype=float)
    col_da_weight = np.zeros((old_x_n, new_x_n), dtype=float)
    for new_x_index in range(new_x_n):
        # determine linear row-weighting of original pixels to new pixels
        temp_edges = new_x_edges[new_x_index:(new_x_index + 2)]
        pixel_portions = pixel_portion_overlap1D(temp_edges, old_x_edges)
        bin_indices = np.where(pixel_portions > 0.)
        bin_weights = pixel_portions[bin_indices]
        col_da_weight[bin_indices, new_x_index] = bin_weights
        # multiply by pixel widths
        bin_weights = bin_weights * old_x_widths[bin_indices]
        # normalize
        bin_weights = bin_weights/bin_weights.sum()
        # store indices and weights for each column
        column_weight_mat[bin_indices, new_x_index] = bin_weights

    # prepare data for weighted-averaging
    full_data = old_mat
    no_data_index = full_data == no_data_val
    full_data[no_data_index] = 0.
    # apply the row and column reduction by matrix multiplication
    row_reduced_data = np.matmul(row_weight_mat, full_data)
    reduced_data = np.matmul(row_reduced_data, column_weight_mat)

    if da is None:
        # create a matrix of ones that matches dimension of old_mat
        da = np.ones(old_mat.shape)

    # also calculate da in the new grid (local linear approximation)
    reduced_grid_da = np.matmul(np.matmul(row_da_weight, da), col_da_weight)
    no_data_da = da.copy()
    no_data_da[no_data_index] = 0.
    reduced_no_data_da = np.matmul(np.matmul(row_da_weight, no_data_da),
                                   col_da_weight)
    # use the area ratio to improve intensity estimate at data boundaries (and
    # better estimate the boundary)
    da_ratio = reduced_no_data_da/reduced_grid_da
    new_no_data_index = da_ratio < 0.5
    reduced_data[new_no_data_index] = no_data_val
    reduced_data[~new_no_data_index] = reduced_data[~new_no_data_index]/ \
        da_ratio[~new_no_data_index]

    return reduced_data


def pixel_portion_overlap1D(edges, new_edges):
    # function to calc new pixel overlap with a single pixel
    """

    :param edges: list-like with two entries (sorted increasing)
           Original pixel edges
    :param new_edges: list-like (sorted increasing)
           New pixel edges
    :return: np.ndarray vector with length=len(new_edges)-1
             The portion of each new pixel that overlaps the original pixel.
    """
    # initiate results vector
    n_bins = len(new_edges) - 1
    out_vec = np.zeros(n_bins)

    left_edges  = new_edges[0:-1]
    right_edges = new_edges[1:]
    left_in     = left_edges  < edges[1]
    left_out    = left_edges  < edges[0]
    right_in    = right_edges > edges[0]
    right_out   = right_edges > edges[1]

    # for extremely large arrays (approx n_bins>2E5), this is probably faster
    # temp_index = np.searchsorted(left_edges, edges[1], side='right')
    # left_in = np.zeros(n_bins, dtype=bool)
    # left_in[:temp_index] = True

    consume_bin = left_out & right_out
    # check for single new pixel that surrounds old pixel
    if any(consume_bin):
        # calculate portion of surrounding bin that overlaps old pixel
        out_vec[consume_bin] = (edges[1] - edges[0])/(right_edges[consume_bin] -
                                                      left_edges[consume_bin])
        return out_vec

    # check left overlap for partial overlap
    left_overlap = left_out & right_in
    out_vec[left_overlap] = (right_edges[left_overlap] - edges[0]) / \
                            (right_edges[left_overlap] - left_edges[left_overlap])

    # check for partial right overlap
    right_overlap = right_out & left_in
    out_vec[right_overlap] = (edges[1] - left_edges[right_overlap]) / \
                             (right_edges[right_overlap] - left_edges[right_overlap])

    # remaining overlap pixels fall inside original pixel
    full_overlap = ~left_out & ~right_out
    out_vec[full_overlap] = 1.

    return out_vec
