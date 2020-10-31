## Jacobi Constant
# =============================================================================
function JC(r, v, μ₃)
    # jacobiConstant(r,v,μ₃)
    # r is [x;y;z] in nondimensional units
    # v is [x';y';z'] in nondimensional units, with as many columns as you like
    # μ₃ is 3-body mu, or m2/(m1+m2)
    # outputs row vector of C-values
    #
    # Alex Hoffman
    # 10/31/2020
    # notation from AAE632 notes

    x, y, z = r[1:1, :], r[2:2, :], r[3:3, :]
    xd, yd, zd = v[1:1, :], v[2:2, :], v[3:3, :]
    d = sqrt.((x .+ μ₃) .^ 2 .+ y .^ 2 .+ z .^ 2)
    r = sqrt.((x .- 1 .+ μ₃) .^ 2 .+ y .^ 2 .+ z .^ 2)
    C =
        -(xd .^ 2 .+ yd .^ 2 .+ zd .^ 2) .+ x .^ 2 .+ y .^ 2 .+
        2(1 - μ₃) ./ d .+ 2μ₃ ./ r

    return C
end

function xtest(a, b)
    return a
end

## Directions Cosine Matrix
# =============================================================================
function DCM(angle, XYZ, isPassive)
    # rot = DCM(angle,XYZ,isPassive)
    # angle is in radians
    # XYZ is 'x' 'y' or 'z'
    # isPassive is a logical stating whether you want a passive rotation
    #   (which rotates the reference frame, not the vector)
    # output is the direction cosines matrix
    # you should use passive when creating a transformation matrix from euler angles
    #
    # Alex Hoffman
    # 10/30/2020

    if XYZ == 'x'
        rot = [1 0 0; 0 cos(angle) sin(angle); 0 -sin(angle) cos(angle)]
    elseif XYZ == 'y'
        rot = [cos(angle) 0 -sin(angle); 0 1 0; sin(angle) 0 cos(angle)]
    elseif XYZ == 'z'
        rot = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1]
    else
        error("Rotation input XYZ must be 'x','y', or 'z'.")
    end
    if ~isPassive
        return transpose(rot)
    else
        return rot
    end
end

## CR#BP Equations of Motion
# =============================================================================
function CR3BP_EOM(state,μ₃,xyzddot)
    # calculates the accelerations in 3BP
    #
    # Alex Hoffman
    # 10/31/2020
    x, y, z, xd, yd, zd = state
    r = state[1:3]
    if xyzddot == 'x'
        dd = 2yd + U_star(r,μ₃,'x')
    elseif xyddot == 'y'
        dd = -2xd + U_star(r,μ₃,'y')
    elseif xyzddot == 'z'
        dd = U_star(r,μ₃,'z')
    else
        error("xyzddot must be 'x', 'y', or 'z'.")
    end
    return dd
end

## Pseudopotential
# =============================================================================
function U_star(r,μ₃,deriv)
    # return U_star(r,mu3,deriv)
    # solves 3-Body U* pseudo-potential and its derivatives
    # X,Y,Z are the position coordinates
    # μ₃ is 3-body mu
    # deriv is 'U','i','ij', where i and j are either 'x','y','z'
    # this represents evaluation of the function and its partial derivatives
    x,y,z = r
    mu3 = μ₃
    a = length(deriv)
    if deriv[1] == 'U'
        U = 0.5*(x^2 + y^2) + (1-mu3)/sqrt((x+mu3)^2+y^2+z^2) + mu3/sqrt((x-1+mu3)^2+y^2+z^2)
        # turn U into symbolic with x,y,z to get derivatives
        return U
    elseif deriv[1] == 'x'
        if a == 1
            Ux = x + ((2*mu3 + 2*x)*(mu3 - 1))/(2*((mu3 + x)^2 + y^2 + z^2)^(3/2)) - (mu3*(2*mu3 + 2*x - 2))/(2*((mu3 + x - 1)^2 + y^2 + z^2)^(3/2))
            return Ux
        else #second partial
            if deriv[2] == 'x'
                Uxx = (mu3 - 1)/((mu3 + x)^2 + y^2 + z^2)^(3/2) - mu3/((mu3 + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu3*(2*mu3 + 2*x - 2)^2)/(4*((mu3 + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu3 + 2*x)^2*(mu3 - 1))/(4*((mu3 + x)^2 + y^2 + z^2)^(5/2)) + 1
                return Uxx
            elseif deriv[2] == 'y'
                Uxy = (3*mu3*y*(2*mu3 + 2*x - 2))/(2*((mu3 + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu3 + 2*x)*(mu3 - 1))/(2*((mu3 + x)^2 + y^2 + z^2)^(5/2))
                return Uxy
            elseif deriv[2] == 'z'
                Uxz = (3*mu3*z*(2*mu3 + 2*x - 2))/(2*((mu3 + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu3 + 2*x)*(mu3 - 1))/(2*((mu3 + x)^2 + y^2 + z^2)^(5/2))
                return Uxz
            else
                error("incorrect derivative expression")
            end
        end
    elseif deriv[1] == 'y'
        if a == 1
            Uy = y - (mu3*y)/((mu3 + x - 1)^2 + y^2 + z^2)^(3/2) + (y*(mu3 - 1))/((mu3 + x)^2 + y^2 + z^2)^(3/2)
            return Uy
        else #second partial
            if deriv[2] == 'x'
                Uyx = (3*mu3*y*(2*mu3 + 2*x - 2))/(2*((mu3 + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu3 + 2*x)*(mu3 - 1))/(2*((mu3 + x)^2 + y^2 + z^2)^(5/2))
                return Uyx
            elseif deriv[2] == 'y'
                Uyy = (mu3 - 1)/((mu3 + x)^2 + y^2 + z^2)^(3/2) - mu3/((mu3 + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu3 - 1))/((mu3 + x)^2 + y^2 + z^2)^(5/2) + (3*mu3*y^2)/((mu3 + x - 1)^2 + y^2 + z^2)^(5/2) + 1
                return Uyy
            elseif deriv[2] == 'z'
                Uyz = (3*mu3*y*z)/((mu3 + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu3 - 1))/((mu3 + x)^2 + y^2 + z^2)^(5/2)
                return Uyz
            else
                error("incorrect derivative expression")
            end
        end
    elseif deriv[1] == 'z'
        if a == 1
            Uz = (z*(mu3 - 1))/((mu3 + x)^2 + y^2 + z^2)^(3/2) - (mu3*z)/((mu3 + x - 1)^2 + y^2 + z^2)^(3/2)
            return Uz
        else #second partial
            if deriv[2] == 'x'
                Uzx = (3*mu3*z*(2*mu3 + 2*x - 2))/(2*((mu3 + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu3 + 2*x)*(mu3 - 1))/(2*((mu3 + x)^2 + y^2 + z^2)^(5/2))
                return Uzx
            elseif deriv[2] == 'y'
                Uzy = (3*mu3*y*z)/((mu3 + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu3 - 1))/((mu3 + x)^2 + y^2 + z^2)^(5/2)
                return Uzy
            elseif deriv[2] == 'z'
                Uzz = (mu3 - 1)/((mu3 + x)^2 + y^2 + z^2)^(3/2) - mu3/((mu3 + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu3 - 1))/((mu3 + x)^2 + y^2 + z^2)^(5/2) + (3*mu3*z^2)/((mu3 + x - 1)^2 + y^2 + z^2)^(5/2)
                return Uzz
            else
                error("incorrect derivative expression")
            end
        end
    else
        error("incorrect derivative expression")
    end
end
