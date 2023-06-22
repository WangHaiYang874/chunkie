classdef level

    properties(Dependent,Access=public)
        nch
    end

    properties(SetAccess=private)
        r
        d
        h
        adj
        p
        pw
        glw
        star
    end

    methods
        function obj = level(r,d,h,adj,p,pw,glw)
            obj.r = r;
            obj.d = d;
            obj.h = h;
            obj.adj = adj;
            obj.p = p;
            obj.pw = pw;
            obj.glw = glw;
            obj.star = false(size(r,3));
        end

        function obj = buildstar(obj)
            obj.star = obj.flagnear();
        end

        function nch = get.nch(obj)
            nch = size(obj.r,3);
        end

        function arclen = arclength(obj)
            wts = reshape(sqrt(sum(obj.d.^2,1)),16,obj.nch);
            wts = wts.*(obj.glw.*(obj.h'));
            arclen = sum(wts,1);
        end

        function near = flagnear(obj)
            arclen = obj.arclength();
            fac = 1/3;
            arclen = fac*arclen;
            lmax = max(arclen)*2.5;
            near = false(obj.nch,1);

            T = hypoct_uni(obj.r(:,:),lmax);

            for i=1:length(T.nodes)

                % getting indices in the tree
                isrc = T.nodes(i).xi;
                inbor = [T.nodes(T.nodes(i).nbor).xi isrc];

                % transfer the point indices into the panel indices
                isrcpanel = idivide(int32(isrc-1),int32(16)) + 1;
                inborpanel = idivide(int32(inbor-1),int32(16)) + 1;

                for j = 1:length(isrc)
                    
                    srcpanel = isrcpanel(j);
                    
                    if near(srcpanel), continue, end

                    distmin = arclen(srcpanel);

                    for k = 1:length(inbor)
                        targpanel = inborpanel(k);
                        if srcpanel == targpanel, continue ,end
                        if srcpanel == obj.adj(1,targpanel), continue, end
                        if srcpanel == obj.adj(2,targpanel), continue, end

                        rj = obj.r(:,isrc(j));
                        rk = obj.r(:,inbor(k));
                        dist = norm(rj-rk);
                        if dist < distmin
                            near(srcpanel) = true;
                        end
                    end
                end
            end
        end

        function nextlvl = next(obj)
           
            star = find(obj.star);
            nstar = sum(obj.star);

            nextr = zeros(2,16,2*nstar);
            nextd = zeros(2,16,2*nstar);
            nexth = zeros(2*nstar,1);

            nextadj = zeros(2,2*nstar);
            p = obj.p;
            pw = obj.pw;
            glw = obj.glw;

            for i=1:nstar

                istar = star(i);

                ri = obj.r(:,:,istar);
                di = obj.d(:,:,istar);
                hi = obj.h(istar);

                nexth(2*i-1) = hi/2;
                nexth(2*i) = hi/2;

                rx = p * ri(1,:)';
                ry = p * ri(2,:)';
                dx = p * di(1,:)';
                dy = p * di(2,:)';
                
                nextr(1,:,2*i-1) = rx(1:16);
                nextr(1,:,2*i) = rx(17:32);
                nextr(2,:,2*i-1) = ry(1:16);
                nextr(2,:,2*i) = ry(17:32);
                
                nextd(1,:,2*i-1) = dx(1:16);
                nextd(1,:,2*i) = dx(17:32);
                nextd(2,:,2*i-1) = dy(1:16);
                nextd(2,:,2*i) = dy(17:32);

                nextadj(2,2*i-1) = 2*i;
                nextadj(1,2*i) = 2*i-1;
            end

            % inheret other adjacency data
            for i=1:nstar
                istar = star(i);
                iprev = obj.adj(1,istar);
                inext = obj.adj(2,istar);
                if iprev ~= 0
                    if obj.star(iprev)
                        nextadj(1,2*i-1) = 2*find(star == iprev);
                    end
                end
                if inext ~= 0 
                    if obj.star(inext)
                        nextadj(2,2*i) = 2*find(star == inext)-1;
                    end
                end
            end

            % build the new level
            nextlvl = chnk.rcip.level(nextr,nextd,nexth,nextadj,p,pw,glw);
            nextlvl = nextlvl.buildstar();
            
        end


        % function kstar_compressed = compress_k_star(obj,kstar)
        %     P = obj.p;
        %     PWT = obj.pw';
        %     nstar = sum(obj.star);
        %     kstar_compressed = zeros(16*nstar,16*nstar);
        %     for i = 1:nstar
        %     for j=1:nstar
        %         i1 = (i-1)*16+1:i*16;
        %         j1 = (j-1)*16+1:j*16;
        %         i2 = 2*(i-1)*16+1:2*i*16;
        %         j2 = 2*(j-1)*16+1:2*j*16;
        %         kstar_compressed(i1,j1) = PWT * kstar(i2,j2) * P;
        %     end
        %     end
        % end

        % function kstar_inv_compressed = compress_k_star_inv(obj,kstar)
        %     P = obj.p;
        %     PWT = obj.pw';
        %     nstar = sum(obj.star);
        %     kstar_compressed = zeros(16*nstar,16*nstar);
        %     for i = 1:nstar
        %     for j=1:nstar
        %         i1 = (i-1)*16+1:i*16;
        %         j1 = (j-1)*16+1:j*16;
        %         i2 = 2*(i-1)*16+1:2*i*16;
        %         j2 = 2*(j-1)*16+1:2*j*16;
        %         kstar_compressed(i1,j1) = PWT * kstar(i2,j2) * P;
        %     end
        %     end
        % end
    end
end 
