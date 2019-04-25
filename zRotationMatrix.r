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

      R <- r_theta %*% r_omega%*% r_phi
      # Order matters, 3rd angle is projection plane rotation.
    }


    rotation_space <- manip_space %*% R
    colnames(rotation_space) <- colnames(manip_space)
    rownames(rotation_space) <- rownames(manip_space)
    rotation_space
  }

}

#runif(1, max = 2 * pi)
#my_t <- 1
rad_theta <- atan(f_bas[f_mvar, 2]/f_bas[f_mvar, 1])
my_p <- 2
f_msp <- create_manip_space(basis = f_bas, manip_var = f_mvar)
b <- rotate_manip_space(f_msp, theta = my_t, phi = my_p) #omega is Null)
b0 <- rotate_manip_space(f_msp, theta = my_p, phi = my_p, omega = 0)
# comparing 2 Euler angle rotaion with 3 angle.
view_basis(b)
view_basis(o0) # far cry from equal, values too small, orthonormality broken?
is_orthonormal(b) # neither are orthonormal
is_orthonormal(o0)
# changing omega; fixed theta, phi. manip_var = 4. rotating whole plane.
o05 <- rotate_manip_space(f_msp, theta = rad_theta, phi = my_p, omega = .5)
o1  <- rotate_manip_space(f_msp, theta = rad_theta, phi = my_p, omega = 1)
o15 <- rotate_manip_space(f_msp, theta = rad_theta, phi = my_p, omega = 1.5)
o2  <- rotate_manip_space(f_msp, theta = rad_theta, phi = my_p, omega = 2)
o25 <- rotate_manip_space(f_msp, theta = rad_theta, phi = my_p, omega = 2.5)
o3  <- rotate_manip_space(f_msp, theta = rad_theta, phi = my_p, omega = 3)
view_basis(o0)
view_basis(o05)
view_basis(o1)
view_basis(o15)
view_basis(o2)
view_basis(o25)
view_basis(o3)


# Fix omega and theta, change phi. manip_var = 4
my_o <- 0
p1  <- rotate_manip_space(f_msp, theta = rad_theta, phi = 1, omega = my_o)
p15 <- rotate_manip_space(f_msp, theta = rad_theta, phi = 1.5, omega = my_o)
p2  <- rotate_manip_space(f_msp, theta = rad_theta, phi = 2, omega = my_o)
p25 <- rotate_manip_space(f_msp, theta = rad_theta, phi = 2.5, omega = my_o)
p3  <- rotate_manip_space(f_msp, theta = rad_theta, phi = 3, omega = my_o)
view_basis(p1)
view_basis(p15)
view_basis(p2)
view_basis(p25)
view_basis(p3)
# swapping the order of phi and omega, to match order of the variables, not expecting any better.


# Fix phi and omega, change theta. manip_var = 4
my_o <- 0
my_p <- 2
t1  <- rotate_manip_space(f_msp, phi = my_p, theta = 1, omega = my_o)
t15 <- rotate_manip_space(f_msp, phi = my_p, theta = 1.5, omega = my_o)
t2  <- rotate_manip_space(f_msp, phi = my_p, theta = 2, omega = my_o)
t25 <- rotate_manip_space(f_msp, phi = my_p, theta = 2.5, omega = my_o)
t3  <- rotate_manip_space(f_msp, phi = my_p, theta = 3, omega = my_o)
view_basis(t1)
view_basis(t15)
view_basis(t2)
view_basis(t25)
view_basis(t3)
# swapping the order of phi and omega, to match order of the variables, not expecting any better.
