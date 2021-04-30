function sur(xP,X0)
% xP: Density matrix; x0: Density threshold

clf

% [nely,nelx,nelz]=size(xP);

xP(:,[1,end],:)=X0;
xP([1,end],:,:)=X0;
xP(:,:,[1,end])=X0;


%
% xP1=xP;
% for i=1:nely
%     for j=1:nelx
%         for k=1:nelz
%             xP(k,j,nelz+1-i)=xP1(i,j,k);
%         end
%     end
% end

% xP1=xP;
% for i=1:nely
%     for j=1:nelx
%         for k=1:nelz
%             xP(k,j,nelz+1-i)=xP1(i,j,k);
%         end
%     end
% end

% xP = smooth3(xP, 'box', 1);
p = patch(isosurface(xP,X0));
isonormals(xP,p);
%  p.FaceColor = 'red'; % Face color
p.FaceColor = '#007FFF'; % Face color
p.EdgeColor = 'none'; % Edge color

daspect([1 1 1])
axis equal; axis tight; axis off;
set(gcf,'color','w');
% material shiny

%%
% Isometric view
camlight(-100,30,'infinite')

% Left
% camlight(100,30,'infinite')

% Right
% camlight(10,0,'infinite')

lighting gouraud

% lighting flat
%  view([30,30]);
% view([45,45,45]);

% Isometric view
view([-45,130,-45]);

% Left
% view([0,45,0]);

% Right
% view([0,0,45]);

% Top 
% view([-100, 0, 0]);

pause(1e-6);


end



