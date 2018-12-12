library(ggplot2)
library(R6)

GameOfLife = R6Class("GameOfLife",
  public = list(
    ## Initialize function
    initialize = function(n_row, n_col) {
      stopifnot(is.numeric(n_row), length(n_row) == 1);
      stopifnot(is.numeric(n_col), length(n_col) == 1);
      
      private$n_row = n_row;
      private$n_col = n_col;
      private$board_mat =
        matrix(
          0L, nrow = n_row, ncol = n_col,
          dimnames = list(seq_len(n_row), seq_len(n_col))
        );
    },
    ## Override default 'print' behavior
    print = function(...) {
      cat("GameOfLife: \n");
      cat("  Board height: ", private$n_row, "\n", sep = "");
      cat("  Board width: ", private$n_col, "\n", sep = "");
      
      return(invisible(self));
    },
    set_init_state = function(pattern = "random") {
      if (pattern == "random") {
        private$board_mat[] = rbinom(private$n_row * private$n_col, 1, 0.1);
      }
      
      return(invisible(self));
    },
    update_game_board = function() {
      n_neighbor_mat = private$obtain_neighbors_amount(private$board_mat);
      private$board_mat = private$apply_survival_rule(private$board_mat, n_neighbor_mat);
      
      return(invisible(self));
    },
    plot_game_board = function() {
      
      board_df = private$melt_matrix_to_dataframe(private$board_mat);
      
      ggplot(data = board_df) +
        geom_rect(
          mapping = aes(xmin = x - 1, xmax = x, ymin = y - 1, ymax = y, fill = status)
        );
    },
    
    # To be removed later
    get_init_state = function() {
      
      return(private$board_mat);
    }
  ),
  private = list(
    n_row = NULL,
    n_col = NULL,
    board_mat = NULL,
    # Private functions
    obtain_neighbors_amount = function(status_mat) {
      stopifnot(min(status_mat) >= 0, max(status_mat) <= 1);
      
      n_input_row = nrow(status_mat);
      n_input_col = ncol(status_mat);
      n_ext_row = n_input_row + 2;
      n_ext_col = n_input_col + 2;
      
      extended_neighbors_mat = matrix(0L, nrow = n_ext_row, ncol = n_ext_col);
      
      row_index_vec = 1L + seq_len(n_input_row);
      col_index_vec = 1L + seq_len(n_input_col);
      
      ## count neightbors in all 8 directions
      extended_neighbors_mat[row_index_vec - 1, col_index_vec - 1] =
        extended_neighbors_mat[row_index_vec - 1, col_index_vec - 1] + status_mat;
      extended_neighbors_mat[row_index_vec - 1, col_index_vec] =
        extended_neighbors_mat[row_index_vec - 1, col_index_vec] + status_mat;
      extended_neighbors_mat[row_index_vec - 1, col_index_vec + 1] =
        extended_neighbors_mat[row_index_vec - 1, col_index_vec + 1] + status_mat;
      extended_neighbors_mat[row_index_vec, col_index_vec - 1] =
        extended_neighbors_mat[row_index_vec, col_index_vec - 1] + status_mat;
      extended_neighbors_mat[row_index_vec, col_index_vec + 1] =
        extended_neighbors_mat[row_index_vec, col_index_vec + 1] + status_mat;
      extended_neighbors_mat[row_index_vec + 1, col_index_vec - 1] =
        extended_neighbors_mat[row_index_vec + 1, col_index_vec - 1] + status_mat;
      extended_neighbors_mat[row_index_vec + 1, col_index_vec] =
        extended_neighbors_mat[row_index_vec + 1, col_index_vec] + status_mat;
      extended_neighbors_mat[row_index_vec + 1, col_index_vec + 1] =
        extended_neighbors_mat[row_index_vec + 1, col_index_vec + 1] + status_mat;
      
      ## compress extended matrix
      ## borders without corners
      min_row_ind = min(row_index_vec);
      max_row_ind = max(row_index_vec);
      min_col_ind = min(col_index_vec);
      max_col_ind = max(col_index_vec);
      extended_neighbors_mat[max_row_ind, col_index_vec] =
        extended_neighbors_mat[max_row_ind, col_index_vec] + extended_neighbors_mat[1, col_index_vec];
      extended_neighbors_mat[min_row_ind, col_index_vec] =
        extended_neighbors_mat[min_row_ind, col_index_vec] + extended_neighbors_mat[n_ext_row, col_index_vec];
      extended_neighbors_mat[row_index_vec, max_col_ind] =
        extended_neighbors_mat[row_index_vec, max_col_ind] + extended_neighbors_mat[row_index_vec, 1];
      extended_neighbors_mat[row_index_vec, min_col_ind] =
        extended_neighbors_mat[row_index_vec, min_col_ind] + extended_neighbors_mat[row_index_vec, n_ext_col];
      ## corners
      lower_dim = min(n_input_row, n_input_col);
      inv_lower_dim_row = n_input_row - lower_dim + 1;
      inv_lower_dim_col = n_input_col - lower_dim + 1;
      ext_row_ind = row_index_vec[lower_dim];
      inv_ext_row_ind = row_index_vec[inv_lower_dim_row];
      ext_col_ind = col_index_vec[lower_dim];
      inv_ext_col_ind = col_index_vec[inv_lower_dim_col];
      
      extended_neighbors_mat[ext_row_ind, ext_col_ind] =
        extended_neighbors_mat[ext_row_ind, ext_col_ind] + extended_neighbors_mat[1, 1];
      extended_neighbors_mat[inv_ext_row_ind, ext_col_ind] =
        extended_neighbors_mat[inv_ext_row_ind, ext_col_ind] + extended_neighbors_mat[n_ext_row, 1];
      extended_neighbors_mat[ext_row_ind, inv_ext_col_ind] =
        extended_neighbors_mat[ext_row_ind, inv_ext_col_ind] + extended_neighbors_mat[1, n_ext_col];
      extended_neighbors_mat[inv_ext_row_ind, inv_ext_col_ind] =
        extended_neighbors_mat[inv_ext_row_ind, inv_ext_col_ind] + extended_neighbors_mat[n_ext_row, n_ext_col];
      
      return(extended_neighbors_mat[row_index_vec, col_index_vec]);
    },
    apply_survival_rule = function(status_mat, neighbor_mat) {
      
      is_alive_mat = status_mat == 1;
      is_two_neighbors_mat = neighbor_mat == 2;
      is_three_neighbors_mat = neighbor_mat == 3;
      
      next_status_mat = is_three_neighbors_mat + is_alive_mat * is_two_neighbors_mat;
      
      return(next_status_mat);
    },
    melt_matrix_to_dataframe = function(input_mat) {
      
      n_row = nrow(input_mat);
      n_col = ncol(input_mat);
      
      index_mat = cbind(rep(seq_len(n_row), n_col), rep(seq_len(n_col), each = n_row));
      melted_mat = cbind(index_mat, input_mat[index_mat]);
      colnames(melted_mat) = c("x", "y", "status");
      melted_df = as.data.frame(melted_mat);
      
      return(melted_df);
    }
  )
);

new_game = GameOfLife$new(n_row = 8, n_col = 10);
new_game$print();
new_game$set_init_state();
new_game$plot_game_board();
new_game$update_game_board();

xxx = new_game$get_init_state();
xxx


