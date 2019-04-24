## scratchpad for changign rotation matrix
# https://en.wikipedia.org/wiki/3D_rotation_group
# https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations

library(spinifex)


if (F) #from rotate_manip_space
{
  rotate_manip_space

  theta = 1
  phi   = 2

  s_theta <- sin(theta)
  c_theta <- cos(theta)
  s_phi <- sin(phi)
  c_phi <- cos(phi)
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
}


# https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations

p=3
#m <- diag(p)
a <- 1/1*pi
b <- 1/2*pi
c <- 0 #1/3*pi # 0 is the null hypothesis for current state.
cos_a <- cos(a)
sin_a <- sin(a)
cos_b <- cos(b)
sin_b <- sin(b)
cos_c <- cos(c)
sin_c <- sin(c)

# cos -sin
# sin  cos
r_a <- matrix(c(1,
                0,
                0,  # 3
                0,
                cos_a,
                -sin_a, # 6
                0,
                sin_a,
                cos_a), # 9
              nrow = 3, ncol = 3,
              byrow = TRUE)

r_b <- matrix(c(cos_b,
                0,
                sin_b, # 3
                0,
                1,
                0, # 6
                -sin_a,
                0,
                cos_a), # 9
              nrow = 3, ncol = 3,
              byrow = TRUE)

r_c <- matrix(c(cos_c,
                -sin_c,
                0, # 3
                sin_c,
                cos_c,
                0, # 6
                0,
                0,
                1), # 9
              nrow = 3, ncol = 3,
              byrow = TRUE)

R <- r_a * r_b * r_c