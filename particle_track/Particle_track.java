/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.lang.Object;
import ucar.nc2.dataset.*;
import ucar.nc2.NCdumpW;
import ucar.nc2.Attribute;
//import ucar.nc2.Netcdf;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;
import ucar.ma2.*;

import java.io.IOException;    

/**
 *
 * @author tomdude
 */
public class Particle_track {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
    
        System.out.println("Starting particle tracking program\n");

        // read in velocity field data
        // directory for file to read
        //String dir="fvcom_output\110817_1306_943node_tideriv_its";

        //System.out.println("Running get_grid2\n");
        
        // get_grid2 reads in data from current_0001.nc to have node and element
        // data. Maybe change java version so that it returns a "mesh" structure?
        // OR - just create the centroids directly by calling a function here that 
        // reads from the netcdf file
        //Mesh mesh = get_grid2(dir);
        //System.out.println("Running plot_usmean_onz\n");
        //plot_usmean_onz(dir);

        
        System.out.println("Reading in data\n");
        //String file= "Z:\\"+dir+"\\output\\netcdf\\current_0001.nc";

        /**
         * files are 
         * meshdata.nc: nodexy, uvnode, trinodes, neigh
         * uv_[3-31].nc: u, v
         */ 
        
        String filename = "D:/sharedfolder/May_2010/meshdata.nc";
        NetcdfFile ncfile = null;
        try {
            ncfile = NetcdfFile.open(filename);
            //process( ncfile);
        } catch (IOException ioe) {
            System.err.println("trying to open " + filename);
        } finally { 
            if (null != ncfile) try {
                ncfile.close();
            } catch (IOException ioe) {
                //log("trying to close " + filename, ioe);
            }
        }
        String varName = "nodeXY";
        Variable nodexy = ncfile.findVariable("nodexy");
        Array data = null;
        if (null == nodexy) return;
        try {
            data = nodexy.read();       
            NCdumpW.printArray(data);
        } catch (IOException ioe) {
            System.err.println("trying to read " + varName);
        }
//            catch (InvalidRangeException e) {
//            System.err.println("invalid Range for " + varName);
//        }
	    //assert(lat.getRank() == 1);	// make sure it's 1-dimensional
        int[] varShape = nodexy.getShape();
        double[][] nodeXY = new double[varShape[0]][varShape[1]];	// where to put them
        Index index = data.getIndex();
        for (int i = 0; i < varShape[0]; i++)
        {
            for (int j = 0; j < varShape[1]; j++)
            {
                nodeXY[i][j]=data.getDouble(index.set(i,j));
                System.out.printf("%.3e ",nodeXY[i][j]);
            }
            System.out.printf("\n");
        }
        
        

        
        //	    int[] index = new int[1]; // index array to specify which value
//	    for (int ilat = 0; ilat < nlats; ilat++) {
//		index[0] = ilat;
//		lats[ilat] = lat.getDouble(index);
//		debug("lats[" + ilat + "]: " + lats[ilat]);
//	    }
//	    /* Read units attribute of lat variable */
//	    String latUnits = lat.getAttribute("units").getStringValue();
//	    debug("attribute lats:units: " + latUnits);

        
        //int M=length(nc_varget(file,"x")); //number of nodes
        //int N=length(nc_varget(file,"partition")); //number of triangles
        
        
        
        
        // will this work directly? check netcdf file variables
        // TODO make function nc_varget -- same data as mesh.uvnodes
        // TODO actually need to interpolate using x,y and nv (list of corner nodes for each element)
//        int nv[] = new int[N];
//        // TODO the element centroids
//        double centroids[][] = makeCentroidXY(file);
//        
//        // the node coordinates 
//        double nodex[]=nc_varget(file,"x");
//        double nodey[]=nc_varget(file,"y");
//
//        // TODO u and v are 3d in netcdf file - read in correctly
//        double u=nc_varget(file,"u");
//        double v=nc_varget(file,"v");
//
//        // setup farms and boundary
//        int tstart=1000;
//        int tmax=3000;
//        int nparts=500;
//        double dt=621;
//
//        double startlocs[][] = setupStartLocations();
//
//        double islandx[] = new double[6];
//        double islandy[] = new double[6];
//        islandx[0] = nodex[655]; islandy[0] = nodey[655];
//        islandx[1] = nodex[642]; islandy[1] = nodey[642];
//        islandx[2] = nodex[641]; islandy[2] = nodey[641];
//        islandx[3] = nodex[629]; islandy[3] = nodey[629];
//        islandx[4] = nodex[640]; islandy[4] = nodey[640];
//        islandx[5] = nodex[656]; islandy[5] = nodey[656];
//        Polygon island = new Polygon(islandx,islandy,6);
//
//        // make a list of boundary nodes ("outline" of the loch)
//        // TODO write setup_bnodes in java
//        int bnodes[]=list_bnodes();
//        double bnodesx[] = new double[bnodes.length];
//        double bnodesy[] = new double[bnodes.length];
//        Polygon boundary = new Polygon(bnodesx,bnodesy,bnodesx.length);
//
//        // boundary.contains(......)
//
//
//
//
//
//        int startcells[] = new int[startlocs[0].length];
//        for (int i = 0; i < startcells.length; i++)
//        {
//            startcells[i]=0;
//            // TODO rewrite whichelement in java - does this need to allow a list of elements to be read in?
//            startcells[i]=whichelement(startlocs[0][i],startlocs[1][i],
//                    nodex,nodey,centroids);
//            // TODO check whether an in polygon method is available in java
//            int inmesh=inpolygon(startlocs[0][i],startlocs[1][i],bnodes);
//            int inisland=inpolygon(startlocs[0][i],startlocs[1][i],island);
//            if (startcells[i]==0)
//            {
//                System.out.printf("INITIAL POINT %d OUTSIDE MESH\n",i);
//            }
//            //System.out.println("//d //d //d ",startcells(i),inmesh,inisland);
//        }
//
//        // make an array of neighbouring elements to each element
//        // TODO rewrite neighbourcells
//        //neighbours=neighbourcells(centroids);
//
//        // the particle tracking
//        System.out.println("Particle tracking starting....\n");
//        startnum=1;
//
//nparts=length(boundary);
//tot_psteps=nparts*nrecords;
//
//if (startnum==1)
//    % an array to save the number of "particle-timesteps" in each cell
//    psteps=zeros(N,1);
//    % array to save source, destination, and transport time for each particle
//    particle_info=zeros(nparts,3);
//    elemPart=zeros(nparts,1);
//end
//
//% setup particles
//particles(nparts,1) = particle(boundary(length(boundary),2),boundary(length(boundary),3));
//if (1==2)
//    xstart=zeros(length(boundary),1);
//    ystart=zeros(length(boundary),1);
//    startElem=zeros(length(boundary),1);
//    for i=1:nparts
//        %fprintf('PARTICLE %d\n',i);
//
//        %startid=int16(length(boundary)*rand(1));
//        startid=i;
//
//        behaviour=1;
//
//        xstart(i)=boundary(startid,2);
//        ystart(i)=boundary(startid,3);
//        %fprintf('start location %d = %d %.4e %.4e\n',i,startid,xstart,ystart);
//        
//    %     inside=inpolygon(xstart,ystart,bnodes(:,1),bnodes(:,2));
//    %     if (inside==0)
//    %         fprintf('INITIAL POINT OUTSIDE MESH\n')
//    %     end
//
//        %
//        %elemPart=whichelement(xstart,ystart,(1:length(mesh.trinodes)),mesh.nodexy,mesh.trinodes)
//
//        % this boundary location is not actually in the mesh/an element, so set
//        % new particle location to centre of nearest element.
//
//        
//        closest=nearest_centroid(xstart(i),ystart(i),centroids);
//
//        xstart(i)=mesh.uvnode(closest,1);
//        ystart(i)=mesh.uvnode(closest,2);
//        startElem(i)=whichelement(xstart(i),ystart(i),(1:length(mesh.trinodes)),mesh.nodexy,mesh.trinodes);
//    end
//    %elemPart
//end
//
//for i=1:nparts
//    particles(i)=particle(xstart(i),ystart(i));
//    particle_info(i,1)=startid;
//    %particles(i).x=startx(i);
//    %particles(i).y=starty(i);
//    elemPart(i)=startElem(i);
//end
//
//% ------------------- loop 2 = timestep ----------------------------
//%fprintf('Starting time loop\n');
//
//
//count0=0;
//count1=0;
//count2=0;
//count3=0;
//count4=0;
//
//minDistTrav=10000000;
//maxDistTrav=0;
//
//tic
//for day=firstday:lastday
//    fprintf('day %d',day);
//    % clear any old data
//    clear FVCOM1
//    % load the new data file. this puts variables straight into the
//    % workspace
//    if (strcmp(datatype,'mat')==true)
//        filename=fvcomfiles(day,1).name;
//        load(filename);
//        u = FVCOM1.u;
//        v = FVCOM1.v;
//        lasttime=24;
//    end
//    
//    
//    lasttime=24;
//
//    for tt=1:lasttime
//        fprintf('\n%d\n',tt);
//%             if (mod(tt,100)==0) 
//%                 fprintf('%d ',tt); 
//%             end
//
//        dep=setParticleDepth(behaviour,tt);
//        uchunk=squeeze(u(tt,dep,:));
//        vchunk=squeeze(v(tt,dep,:));
//        velocities=[uchunk vchunk];
//
//        for st=1:stepsPerStep
//            %fprintf(',\n');
//            
//
//            for i=1:nparts
//                %fprintf('PARTICLE %d\n',i);
//                if (particles(i).arrived==0)
//
//                    if (veltype==1)
//                        % 1. find neighbouring cells to compute velocity
//                        closest=nearest_centroid(particles(i).x,particles(i).y,centroids);
//                        %closest=whichelement(particles(i).x,particles(i).y,(1:length(grid.mesh2.trinodes)),grid.mesh2.nodexy,grid.mesh2.trinodes);
//                        % 2. calculate velocity applied to particle
//                        water_U=velocities(closest,1);
//                        water_V=velocities(closest,2);
//                        %fprintf('water u = %e v = %e\n',water_U,water_V);
//                    end
//                    if (veltype==2)
//                        [water_U,water_V,closest]=velocity_nearest5centroids(p,centroids,velocities);
//                        %fprintf('water u = %e v = %e\n',water_U,water_V);
//                    end
//                    % 
//                    if (veltype==3)
//                        % look in the current element and its neighbours
//                        %elemPart(i)=whichelement(particles(i).x,particles(i).y,neighbours(elemPart(i),:),mesh.nodexy,mesh.trinodes);
//                        
//                        % this bit combines the finding of element for velocity
//                        % computation with the location update, preventing the
//                        % need to check whether the particle is within the
//                        % region at the end of the timestep.
//                        % BASIC IDEA - leave location (x,y,cell) unchanged if
//                        % the particle has jumped out of the water
//
//                        %fprintf('%d',elemPart(i));
//                        water_U=velocities(elemPart(i),1);
//                        water_V=velocities(elemPart(i),2);
//                        closest=elemPart(i);
//                    end
//                    % save location information "start of timestep"
//                    psteps(elemPart(i))=psteps(elemPart(i))+1;
//
//                    %[water_U,water_V]=calc_vel1(particles(i),[uchunk vchunk]);
//                    % 3. Calculate diffusion    
//                    %rand('twister',sum(100*clock)); %resets it to a different state each time.
//                    diff_X = sqrt(0.01*dt);
//                    diff_Y = sqrt(0.01*dt);    %+/- is random so direction doesn't matter
//                    [behave_U,behave_V]=behaveVelocity(behaviour);
//                    % 4. update particle location
//                    newlocx=particles(i).x+dt*water_U+dt*behave_U+normrnd(0,diff_X); % simplest possible "Euler"
//                    newlocy=particles(i).y+dt*water_V+dt*behave_V+normrnd(0,diff_Y);
//                    
//
//                    %fprintf('+');
//                    whereami=whichelement(newlocx,newlocy,elemPart(i),mesh.nodexy,mesh.trinodes);
//                    count0=count0+1;
//                    if (whereami==0)
//                        checkfirst=[elemPart(i),CVM.NBE(elemPart(i),find(CVM.NBE(elemPart(i),:)))];
//                        count1=count1+1;
//                        whereami=whichelement(newlocx,newlocy,checkfirst,mesh.nodexy,mesh.trinodes);
//                        % if fails, look in nearest 10 (id numerical)
//                        if (whereami==0)
//                            %fprintf('a');
//                            %checkfirst
//                            count2=count2+1;
//                            whereami=whichelement(newlocx,newlocy,max(1,elemPart(i)-5):min(length(mesh.trinodes),elemPart(i)+5),mesh.nodexy,mesh.trinodes);
//                            % if fails, look in nearest 500 (id numerical)
//                            if (whereami==0)
//                                %fprintf('b');
//                                %checkfirst
//                                count3=count3+1;
//                                whereami=whichelement(newlocx,newlocy,max(1,elemPart(i)-250):min(length(mesh.trinodes),elemPart(i)+250),mesh.nodexy,mesh.trinodes);
//                                % if this fails, look in all elements
//                                if (whereami==0)
//                                    %fprintf('c');
//                                    count4=count4+1;
//                                    whereami=whichelement(newlocx,newlocy,(1:length(mesh.trinodes)),mesh.nodexy,mesh.trinodes);
//                                    %elemPart(i)=nearest_centroid(particles(i).x,particles(i).y,centroids);
//                                end
//                            end
//                        end
//                    end
//                    
//                    
//                    if (whereami~=0)
//                        distTrav=sqrt((particles(i).x-newlocx).^2+(particles(i).y-newlocy).^2);
//                        if (distTrav>maxDistTrav)
//                            maxDistTrav=distTrav;
//                        end
//                        if (distTrav<minDistTrav)
//                            minDistTrav=distTrav;
//                        end
//                        particles(i).x=newlocx;
//                        particles(i).y=newlocy;
//                        elemPart(i)=whereami;
//                    end
//                    if (whereami==0)
//                        closest=nearest_centroid(particles(i).x,particles(i).y,centroids);
//                        %fprintf('x%d',closest);
//                        particles(i).x=mesh.uvnode(closest,1);
//                        particles(i).y=mesh.uvnode(closest,2);
//                        elemPart(i)=closest;
//                    end
//                
//                    % set particle to become able to settle after a predefined time
//                    if ((day-1)*24+tt>viabletime)
//                        particles(i).viable=1;
//                    end
//                    // if able to settle, is it close to a possible settlement
//                    // location?
//                    if (particles(i).viable==1)
//                        for bnode=1:length(boundary)
//                            dist=sqrt((particles(i).x-boundary(bnode,2)).^2+(particles(i).y-boundary(bnode,3)).^2);
//                            if (dist < 250)
//                                particles(i).arrived=1;
//                                particle_info(i,2)=bnode;
//                                particle_info(i,3)=((day-1)*24+tt);
//                                break
//                            end
//                        end    
//                    end
//                    % break time loop
//    %                 if (particles(i).arrived==1)
//    %                     break;
//    %                 end
//                end
//            end
//            % end of particle loop
//        end
//    end
//    % break day loop
//%     if (particles(i).arrived==1)
//%         break;
//%     end
//end   
//fprintf('element search counts: %d %d %d %d %d\n',count0,count1,count2,count3,count4);
//fprintf('transport distances: min = %.4e, max = %.4e', minDistTrav, maxDistTrav);
//
//       
//        //profile viewer
//        //p = profile("info");
//        //profsave(p,"profile_results")
//        // scale psteps
//        tot_psteps=tmax*nparts;
//        for j=1:length(startlocs)
//            for i=1:N
//                psteps2(i,j)=psteps(i,j)/tot_psteps;
//            end
//        end
//       
    }

    public double[][] setupStartLocations()
    {
        int nsites = 9;
        double startlocs[][]= new double[2][nsites];
        
        startlocs[0][0]=357420; startlocs[1][0]=6217200; 
        startlocs[0][1]=361834; startlocs[1][1]=6223063;
        startlocs[0][2]=353078; startlocs[1][2]=6206339;
        startlocs[0][3]=354246; startlocs[1][3]=6194759; 
        startlocs[0][4]=352745; startlocs[1][4]=6201735;
        startlocs[0][5]=348880; startlocs[1][5]=6199380;
        startlocs[0][6]=354969; startlocs[1][6]=6193169;
        startlocs[0][7]=348606; startlocs[1][7]=6204475;
        startlocs[0][8]=352401; startlocs[1][8]=6190933;
        // non-fishfarms
//        double startlocs[][]=[354500 6188000; 
//            355200 6192000;
//            350000 6197000;
//            352800 6196000; 
//            348000 6209000;
//            354000 6204000;
//            357000 6213000;
//            360000 6219000;
//            370000 6227000];
        // "hypothetical" fishfarm used for SAMS newsletter
        //fishfarms(1,1)=351000;
        //fishfarms(1,2)=6195000;
        //"2" is Minard
        //"7" is Portavadie
        return startlocs;
    }

    public double[][] makeCentroidXY()
    {
        double[][] out=new double[1][1];
        return out;
    }

    public void setupOutput()
    {

    }

    public void writeOutput()
    {
    
    }
    
    public double nc_varget(String filename, String variable)
    {
        double out=0;
        return out;
    }

}
