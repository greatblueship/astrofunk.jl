function DCM(angle,XYZ,isPassive)
    # rot = DCM(angle,XYZ,isPassive)
    # angle is in radians
    # XYZ is 'x' 'y' or 'z'
    # isPassive is a logical stating whether you want a passive rotation
    #   (which rotates the reference frame, not the vector)
    # output is the direction cosines matrix
    # you should use passive when creating a transformation matrix from euler angles

    if XYZ == 'x'
        rot = [1 0 0; 0 cos(angle) sin(angle); 0 -sin(angle) cos(angle)]
    elseif XYZ == 'y'
        rot =  [cos(angle) 0 -sin(angle); 0 1 0; sin(angle) 0 cos(angle)]
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
