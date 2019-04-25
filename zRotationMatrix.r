## scratchpad for changign rotation matrix
# https://en.wikipedia.org/wiki/3D_rotation_group
# https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations

library(spinifex)

set.seed(20190425)
f_dat  <- tourr::rescale(flea[,1:6])
f_cat  <- factor(flea$species)
f_bas  <- basis_random(6,2)
f_mvar <- 4

view_manip_space(basis = f_bas,manip_var = f_mvar)
#play_manual_tour(data = f_dat, basis = f_bas, manip_var = f_mvar,
#                 col = col_of(flea$species), axes = "bottomleft")

load_funcs <- T
if (load_funcs) # overwrites the spinifex functions: rotate_manip_space and manual_tour
{
  #rotate_manip_space
  rotate_manip_space <- function(manip_space,
                                    theta=NULL,
                                    phi=NULL,
                                    omega=NULL # if null, use Cook's 2 angle rotation matrix.
  )
  {
    # if(is.null(theta)) theta = 1
    # if(is.null(phi))   phi   = 2
    s_theta <- sin(theta)
    c_theta <- cos(theta)
    s_phi <- sin(phi)
    c_phi <- cos(phi)
    if(!is.null(omega))
      {s_omega <- sin(omega)
      c_omega <- cos(omega)}

    if(is.null(omega)){ # if omega not passed, use Cook defined (2 angle) rotation.
      R <- matrix(c(c_theta^2 * c_phi + s_theta^2,
                    -c_theta * s_theta * (1 - c_phi),
                    -c_theta * s_phi,  # 3
                    -c_theta * s_theta * (1 - c_phi),
                    s_theta^2 * c_phi + c_theta^2,
                    -s_theta * s_phi, # 6
                    c_theta * s_phi,
                    s_theta * s_phi, c_phi), # 9
                  nrow = 3, ncol = 3,
                  byrow = TRUE)
    } else { #if omega is passed, use 3 euler angle rotation
      r_theta <- matrix(c(1,
                      0,
                      0,  # 3
                      0,
                      c_theta,
                      -s_theta, # 6
                      0,
                      s_theta,
                      c_theta), # 9
                    nrow = 3, ncol = 3,
                    byrow = TRUE)

      r_phi <- matrix(c(c_phi,
                      0,
                      s_phi, # 3
                      0,
                      1,
                      0, # 6
                      -s_phi,
                      0,
                      c_phi), # 9
                    nrow = 3, ncol = 3,
                    byrow = TRUE)

      r_omega <- matrix(c(c_omega,
                      -s_omega,
                      0, # 3
                      s_omega,
                      c_omega,
                      0, # 6
                      0,
                      0,
                      1), # 9
                    nrow = 3, ncol = 3,
                    byrow = TRUE)

      R <- r_theta * r_phi * r_omega
    }


    rotation_space <- manip_space %*% R
    colnames(rotation_space) <- colnames(manip_space)
    rownames(rotation_space) <- rownames(manip_space)
    rotation_space
  }

  # # paste over spinifex function
  # manual_tour <- function (basis = NULL, manip_var, theta = NULL, phi_min = 0,
  #                          phi_max = 0.5 * pi, omega=NULL, angle = 0.05)
  # {
  #   basis <- as.matrix(basis)
  #   p <- nrow(basis)
  #   d <- ncol(basis)
  #   manip_space <- create_manip_space(basis = basis, manip_var = manip_var)
  #   if (is.null(theta))
  #     theta <- atan(basis[manip_var, 2]/basis[manip_var, 1])
  #   phi_start <- acos(sqrt(basis[manip_var, 1]^2 + basis[manip_var, 2]^2))
  #   stopifnot(phi_min <= phi_start & phi_max >= phi_start)
  #   find_path <- function(start, end) {
  #     mvar_xsign <- -sign(basis[manip_var, 1])
  #     start <- mvar_xsign * (start - phi_start)
  #     end <- mvar_xsign * (end - phi_start)
  #     dist <- abs(end - start)
  #     steps <- round(dist/angle)
  #     angle <- dist/steps
  #     sign <- ifelse(end > start, 1, -1)
  #     seq(from = start, to = end, by = sign * angle)
  #   }
  #   phi_path <- c(find_path(start = phi_start, end = phi_min),
  #                 find_path(start = phi_min, end = phi_max), find_path(start = phi_max,
  #                                                                      end = phi_start))
  #   n_frames <- length(phi_path)
  #   basis_set <- array(NA, dim = c(p, d, n_frames))
  #   for (i in 1:n_frames) {
  #     thisFrame <- my_rotate_manip_space(manip_space = manip_space,
  #                                        theta = theta, phi = phi_path[i],
  #                                        omega = omega)
  #     basis_set[, , i] <- thisFrame[, 1:d]
  #   }
  #   attr(basis_set, "manip_var") <- manip_var
  #   basis_set
  # }

}

#runif(1, max = 2 * pi)
my_t <- 1
my_p <- 2
f_msp <- create_manip_space(basis = f_bas, manip_var = f_mvar)
b <- rotate_manip_space(f_msp, theta = my_t, phi = my_p) #omega is Null)
b0 <- rotate_manip_space(f_msp, theta = my_p, phi = my_p, omega = 0)
# comparing 2 Euler angle rotaion with 3 angle.
view_basis(b)
view_basis(b0) # far cry from equal, values too small, orthonormality broken?
is_orthonormal(b) # neither are orthonormal
is_orthonormal(b0)
# changing omega; fixed theta, phi. manip_var = 4
b05 <- rotate_manip_space(f_msp, theta = my_p, phi = my_p, omega = .5)
b1  <- rotate_manip_space(f_msp, theta = my_p, phi = my_p, omega = 1)
b15 <- rotate_manip_space(f_msp, theta = my_p, phi = my_p, omega = 1.5)
b2  <- rotate_manip_space(f_msp, theta = my_p, phi = my_p, omega = 2)
b25 <- rotate_manip_space(f_msp, theta = my_p, phi = my_p, omega = 2.5)
b3  <- rotate_manip_space(f_msp, theta = my_p, phi = my_p, omega = 3)
view_basis(b0)
view_basis(b05)
view_basis(b1)
view_basis(b15)
view_basis(b2)
view_basis(b25)
view_basis(b3)


# Fix omega and theta, change phi.  manip_var = 4
my_t <- 1
my_o <- 4
p1  <- rotate_manip_space(f_msp, theta = my_p, phi = 1, omega = my_o)
p15 <- rotate_manip_space(f_msp, theta = my_p, phi = 1.5, omega = my_o)
p2  <- rotate_manip_space(f_msp, theta = my_p, phi = 2, omega = my_o)
p25 <- rotate_manip_space(f_msp, theta = my_p, phi = 2.5, omega = my_o)
p3  <- rotate_manip_space(f_msp, theta = my_p, phi = 3, omega = my_o)
view_basis(p1)
view_basis(p15)
view_basis(p2)
view_basis(p25)
view_basis(p3)

