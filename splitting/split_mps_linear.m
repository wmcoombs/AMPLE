function [mpData,split_collection] = split_mps_linear(mesh,mpData,directions)
to_split = zeros(length(mpData),1);
nmp = length(mpData);
nmp_new = nmp+1;
split_collection = zeros(2*sum(directions>0),1);
scol = 1;
for mp=1:nmp 
    if mp <= length(directions)
    direction = directions(mp);
    if direction ~= 0
        mpData(nmp_new) = mpData(mp);
        position_disp = mpData(mp).lp(direction) * 0.5;
        new_vp = mpData(mp).vp * 0.5;
        new_vp0 = mpData(mp).vp0 * 0.5;
        new_mpM  = mpData(mp).mpM * 0.5;
        
        ids = [mp,nmp_new];
        normals = [-1,1];
        for i = 1:2
            cid = ids(i);
            split_collection(scol) = cid;
            scol = scol + 1;
            normal = normals(i);
            mpData(cid).mpC(direction) = mpData(cid).mpC(direction) + (position_disp*normal);
            mpData(cid).lp(direction) = mpData(cid).lp(direction)*0.5;
            mpData(cid).lp0(direction)= mpData(cid).lp0(direction)*0.5;
            mpData(cid).vp = new_vp;
            mpData(cid).vp0 = new_vp0;
            mpData(cid).mpM = new_mpM;
            x = mpData(cid).lp(1);y = mpData(cid).lp(2);                                          % size of the material point in x and y
            cmat = mpCorners(x,y) + repmat(mpData(cid).mpC,4,1);                 
            mpData(cid).C      = cmat;

            if ((-x + mpData(cid).mpC(1))<=0)
                x = mpData(cid).mpC(1)-10^-7;
            end
            max_X = mesh.size(1);
            if ((x + mpData(cid).mpC(1))>max_X)
                x = (max_X - mpData(cid).mpC(1))-10^-7;
            end
            mpData(cid).lp(1) = x;
    
    
            if ((-y + mpData(cid).mpC(2))<=0)
                y = mpData(cid).mpC(2)-10^-7;
            end
            max_Y = mesh.size(2);
            if ((y + mpData(cid).mpC(2))>max_Y)
                y = (max_Y - mpData(cid).mpC(2))-10^-7;
            end
            mpData(cid).lp(2) = y;
            mpData(cid).C = mpCorners(x,y) + repmat(mpData(cid).mpC,4,1);
        end
        nmp_new =  nmp_new + 1;
    end
    end
end

