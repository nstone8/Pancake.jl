module Pancake

using Tessen, Delica, Unitful, Ahmes, Statistics, Primes
import Tessen:HatchLine, pointalong, intersections

export createconfig, scaffold

"""
```julia
createconfig([filename])
```
Write an example config file to `filename`. If `filename` is omitted, the file
will be written to `"config.jl"`.

# Configuration Parameters
- lscaf: scaffold length
- wscaf: scaffold width
- dfield: calibrated FOV for the objective being used
- hbottom: distance from the bottom of the post to the bottom of the beams
- htop: distance from the bottom of the beams to the top of the bumpers
- chamferbottom: angle at which the bottoms of the posts and beams should be chamfered
- chamfertop: angle at which the tops of the posts and beams should be chamfered
- cutangle: angles at which blocks should be cut to avoid shadowing
- overlap: amount that neighboring blocks should overlap to ensure they are connected
- nhammock: number of hammock slices to write
- dhammockslice: distance between hammock slices
- dslice: slicing distance
- wbumper: width of the bumpers
- fillet: fillet radius on xz and yz crossections of the posts and bumpers
- wpost: post width
- hbeam: beam height
- lbeammax: maximum total beam length
- maxseglength: maximum length of a single beam segment
- keygap: closest point between the two halves of a beam before they are closed with a keystone
- dhammockhatch: hammock hatching distance
- dhatch: hatching distance for posts, beams and bumpers
- laserpower: laser power for posts, beams and bumpers
- scanspeed: scan speed for posts, beams and bumpers
- stagespeed: max stage speed
- interfacepos: how far into the substrate writing should begin
- hamscanspeed: scan speed for hammocks
- hamlaserpower: laser power for hammocks
"""
function createconfig(filename="config.jl")
    config = """Dict(
        :lscaf => 6u"mm",
        :wscaf => 4u"mm",
        :dfield => 1750u"µm",
        :hbottom => 30u"µm",
        :htop => 50u"µm",
        :chamferbottom => pi/6,
        :chamfertop => pi/6,
        :cutangle => pi/6,
        :overlap => 10u"µm",
        :dslice => 1u"µm",
        :nhammock => 3,
        :dhammockslice => 300u"nm",
        :wbumper => 200u"µm",
        :fillet => 30u"µm",
        :wpost => 25u"µm",
        :hbeam => 10u"µm",
        :lbeammax => 120u"µm",
        :maxseglength => 30u"µm",
        :keygap => 20u"µm",
        :dhammockhatch => 1u"µm",
        :dhatch => 300u"nm",
        :laserpower => 100u"mW",
        :scanspeed => 100000u"µm/s",
        :stagespeed => 100u"µm/s",
        :interfacepos => 10u"µm",
        :hamscanspeed => 10000u"µm/s",
        :hamlaserpower => 100u"mW",
    )
    """
    open(filename,"w") do io
        print(io,config)
    end
end

function bumperedgecoords(;kwargs...)
    #copy out the variables we want
    hbottom = kwargs[:hbottom]
    htop = kwargs[:htop]
    chamfertop = kwargs[:chamfertop]
    chamferbottom = kwargs[:chamferbottom]
    wbumper = kwargs[:wbumper]
    fillet = kwargs[:fillet]
    #figuring out how to fillet everything in crossection seems hard enough that I'm just going
    #to use Tessen hatching machinery to do it for me.
    #this function will return a closure which gives the bottom and top coordinates of a bumper
    #as a function of z position

    #draw our crossection on the xz plane
    verts=[
        [-wbumper/2,hbottom],
        [-wbumper/2 + htop*tan(chamfertop),hbottom+htop],
        [wbumper/2 - htop*tan(chamfertop),hbottom+htop],
        [wbumper/2,hbottom],
        [wbumper/2 - hbottom*tan(chamferbottom),0u"µm"],
        [-wbumper/2 + hbottom*tan(chamferbottom),0u"µm"]
    ]
    nofillet = polycontour(verts)
    filleted = polycontour(verts,fillet)
    function(zcoord)
        #place a hatchline off of the contour to the left and at our desired height
        hl=HatchLine([-wbumper,zcoord],[1,0]u"µm")
        #our first coordinate is the first intersection of hl with nofillet
        firstpara = sort(intersections(nofillet,hl))[1]
        #our second coordinate is the second intersection of hl with filleted
        secondpara = sort(intersections(filleted,hl))[2]
        map([firstpara,secondpara]) do p
            #convert to real world coordinates in µm
            coords = pointalong(hl,p)
            #we're interested in the x coordinate, add the dimensions on
            coords[1] * u"µm"
        end
    end
end

#helper function to build a series of segments to get us from startx to endx
#we will assume startx and endx do not include overlap on either side
function bumpersegtrain(startx,endx,startz,endz,bec;kwargs...)
    #to connect the left and right sides, we will use a bunch of identical segments cut to
    #have overhangs on the left side and placed to overlap with each other
    #get the amount of length that should be devoted to overhangs
    lseg = sqrt(kwargs[:dfield]^2 - kwargs[:wbumper]^2)
    olength = 2*(endz-startz)*tan(kwargs[:cutangle])
    #the length of every slice should be constant, the maximum progress we can make per slice
    #is set by
    lslicemax = lseg - (olength/2) - kwargs[:overlap]
    #the distance we have to cover (including overlapping at either end) is...
    reqdist = endx -  startx + 2*sign(endx-startx)*kwargs[:overlap]

    numsegs = ceil(Int,abs(reqdist/lslicemax))
    distperseg = reqdist/numsegs

    #build a block representing one of these segments with the middle of the bottom of the `source`
    #edge at the local origin

    segslices = map(range(start=startz,step=kwargs[:dslice],stop=endz)) do z
        ycoords = bec(z)
        #we should slant the opposite direction when moving backwards
        sliceoffset = -sign(reqdist)*z*tan(kwargs[:cutangle])
        #draw the edges
        e = [LineEdge([sliceoffset,ycoords[1]],
                      [sliceoffset+distperseg+sign(distperseg)*kwargs[:overlap],ycoords[1]]),
             LineEdge([sliceoffset+distperseg+sign(distperseg)*kwargs[:overlap],ycoords[1]],
                      [sliceoffset+distperseg+sign(distperseg)*kwargs[:overlap],ycoords[2]]),
             LineEdge([sliceoffset+distperseg+sign(distperseg)*kwargs[:overlap],ycoords[2]],
                      [sliceoffset,ycoords[2]]),
             LineEdge([sliceoffset,ycoords[2]],
                      [sliceoffset,ycoords[1]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    #place the first such block overlapping with lbblock
    firstseg = Block(segslices...,origin=[startx - sign(distperseg)*kwargs[:overlap],
                                          0u"µm",startz])
    #we actually want the objective to be centered on the segment during the print, so we need to
    #move the local origin. we will do this by sliding the geometry to the left while using
    #preserveframe, and then slide the resulting object back
    #the -startz translation in the preserveframe step moves the geometry to z=0 in the local
    #frame, the origin of which is at startz in the enclosing frame
    fbsdisp = sign(distperseg)*olength+(distperseg/2)
    blocks = [translate(translate(firstseg,[-fbsdisp,0u"µm",-startz],preserveframe=true),
                        [fbsdisp,0u"µm"])]
    #build the rest by copying
    for _ in 2:numsegs
        push!(blocks,translate(blocks[end],[distperseg,0u"µm"]))
    end
    return blocks
end

function bumper(;kwargs...)
    #The bumper is going to be way longer than `dfield`. The longest segment with width `w`
    #we can write at a time is
    lseg = sqrt(kwargs[:dfield]^2 - kwargs[:wbumper]^2)
    #this includes all of the space we have for overlaps

    #get our width profile as a function of z coordinate
    bec=bumperedgecoords(;kwargs...)
    #we first need to build the bottom of the leftmost end of the bumper. the width of this block will come
    #from `bec`, the left face will be completely rounded off.
    #we would like the right face of the block to taper down to the substrate at `cutangle`
    #we would like this taper to hit the substrate at x=lseg/2, and for the rounded arc of the
    #widest slice (which will have a radius of `wbumper` to be tangent to x=-lseg/2

    #we will build the bumper up to a z height of (hbottom+htop+overlap)/2 on the first pass, and
    #then glue on the top coming back
    zmid=(kwargs[:hbottom]+kwargs[:htop]+kwargs[:overlap])/2
    lbslices = map(range(start=0u"µm",step=kwargs[:dslice],stop=zmid)) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom right corner, this should be lseg/2
        #for the bottom slice
        xbr = (lseg/2) - z*tan(kwargs[:cutangle])
        #get the center of the arc making up the left face
        cleftarc = [-lseg/2 + kwargs[:wbumper]/2, mean(ycoords)]
        #draw the edges
        e = [ArcEdge(cleftarc,wslice/2,pi/2,3pi/2),
             LineEdge([cleftarc[1],ycoords[1]],[xbr,ycoords[1]]),
             LineEdge([xbr,ycoords[1]],[xbr,ycoords[2]]),
             LineEdge([xbr,ycoords[2]],[cleftarc[1],ycoords[2]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    #we want the origin of this block to be at [(-lscaf+lseg)/2,0] so that the whole object is
    #centered on the origin
    lbblock=Block(lbslices...,origin=[(-kwargs[:lscaf]+lseg)/2,0u"µm",0u"µm"])

    #we will now do the bottom of the rightmost block. This will be the mirror image of lbblock
    #but with the cut going the opposite way (i.e. the topmost slice is longest)
    rbslices = map(range(start=0u"µm",step=kwargs[:dslice],stop=zmid)) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom left corner, this should be -lseg/2
        #for the topmost slice
        xbl = (-lseg/2) + (zmid-z)*tan(kwargs[:cutangle])
        crightarc = [lseg/2 - kwargs[:wbumper]/2, mean(ycoords)]
        #draw the edges
        e = [ArcEdge(crightarc,wslice/2,3pi/2,pi/2),
             LineEdge([crightarc[1],ycoords[2]],[xbl,ycoords[2]]),
             LineEdge([xbl,ycoords[2]],[xbl,ycoords[1]]),
             LineEdge([xbl,ycoords[1]],[crightarc[1],ycoords[1]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    #we want the origin of this block at [(lscaf-lseg)/2,0]
    rbblock=Block(rbslices...,origin=[(kwargs[:lscaf]-lseg)/2,0u"µm",0u"µm"])

    #to connect the left and right sides, we will use a bunch of identical segments cut to
    #have overhangs on the left side and placed to overlap with each other
    midbottomblocks = bumpersegtrain(-kwargs[:lscaf]/2 + lseg,
                                     kwargs[:lscaf]/2 - lseg + zmid*tan(kwargs[:cutangle]),
                                     0u"µm",zmid,bec;kwargs...)
    #now we have to build the top. We'll do similar to what we did for the bottom, but we'll
    #make the end pieces sorter so the seams between blocks on the top and bottom are less
    #likely to line up

    zmid2=(kwargs[:hbottom]+kwargs[:htop]-kwargs[:overlap])/2
    rtslices = map(range(start=zmid2,step=kwargs[:dslice],stop=kwargs[:hbottom]+kwargs[:htop])) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom left corner, this should be -lseg/4
        #for the bottom slice
        xbl = (-lseg/4) + (z)*tan(kwargs[:cutangle])
        crightarc2 = [lseg/4 - kwargs[:wbumper]/2, mean(ycoords)]
        #draw the edges
        e = [ArcEdge(crightarc2,wslice/2,3pi/2,pi/2),
             LineEdge([crightarc2[1],ycoords[2]],[xbl,ycoords[2]]),
             LineEdge([xbl,ycoords[2]],[xbl,ycoords[1]]),
             LineEdge([xbl,ycoords[1]],[crightarc2[1],ycoords[1]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    rtblock=Block(rtslices...,origin=[kwargs[:lscaf]/2-lseg/4,0u"µm",zmid2])
    #given that we want the origin of the coordinate system at zmid2, we need to translate
    #the geometry in rtblock to start at z=0 in the local frame
    rtblock=translate(rtblock,[0u"µm",0u"µm",-zmid2],preserveframe=true)

    ltslices = map(range(start=zmid2,step=kwargs[:dslice],stop=kwargs[:hbottom]+kwargs[:htop])) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom right corner, this should be lseg/4
        #for the topmost slice
        xbr = (lseg/4) + (z-kwargs[:htop])*tan(kwargs[:cutangle])
        cleftarc2 = [-lseg/4 + kwargs[:wbumper]/2, mean(ycoords)]
        #draw the edges
        e = [ArcEdge(cleftarc2,wslice/2,pi/2,3pi/2),
             LineEdge([cleftarc2[1],ycoords[1]],[xbr,ycoords[1]]),
             LineEdge([xbr,ycoords[1]],[xbr,ycoords[2]]),
             LineEdge([xbr,ycoords[2]],[cleftarc2[1],ycoords[2]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    ltblock=Block(ltslices...,origin=[-kwargs[:lscaf]/2+lseg/4,0u"µm",zmid2])
    #given that we want the origin of the coordinate system at zmid2, we need to translate
    #the geometry in ltblock to start at z=0 in the local frame
    ltblock=translate(ltblock,[0u"µm",0u"µm",-zmid2],preserveframe=true)

    #build the middle segments in the top
    midtopblocks = bumpersegtrain(kwargs[:lscaf]/2 - lseg/2,
                                  -kwargs[:lscaf]/2 + lseg/2 -
                                      (kwargs[:htop]-zmid2)*tan(kwargs[:cutangle]),
                                  zmid2,kwargs[:hbottom]+kwargs[:htop],bec;kwargs...)
    return SuperBlock(lbblock,midbottomblocks...,rbblock,rtblock,midtopblocks...,ltblock)
end

#build a segmented beam with nsegs segments centered on y=0.
#startx and stopx should be the position of the edge we are bonding to at the z coordinate
#corresponding to the bottom of the beam (i.e. this function will build in overlap)
#gonna go overboard and make this a struct so we can define rotation and translation
struct Beam
    leftsegs
    rightsegs
    keystone
end

function Beam(nsegs::Int,startx,stopx,width;kwargs...)
    #nsegs needs to be odd
    @assert isodd(nsegs)
    #we want the keystone to be as short as possible. Point of closest approach between the
    #two halves is set by `keygap`
    #length of the keystone measured at the beam midline
    lkey = kwargs[:keygap] + kwargs[:hbeam]*tan(kwargs[:cutangle])
    #the first segment needs to be chamfered according to `chamfertop`, we will do the rest at
    #cutangle
    distperseg = (stopx - startx - lkey)/(nsegs-1)
    lseg = distperseg + kwargs[:overlap]
    leftsegs = map(1:((nsegs-1)/2)) do i
        #the first segment needs to be a little longer and chamfered differently
        segpos = startx + distperseg*((2i-1)/2)
        seg = if i==1
            box(lseg + (kwargs[:overlap]/2),width,kwargs[:hbeam],kwargs[:dslice],
                chamfer =[-kwargs[:chamfertop] kwargs[:cutangle]
                          kwargs[:cutangle]    kwargs[:cutangle]])            
        else
            box(lseg,width,kwargs[:hbeam],kwargs[:dslice],
                chamfer =[-kwargs[:cutangle]   kwargs[:cutangle]
                          kwargs[:cutangle]    kwargs[:cutangle]])            
        end
        #seg is currently centered at [0,0,0]. Move it into position (use preserveframe so we don't
        #move the stage
        seg=translate(seg,[segpos,0u"µm",kwargs[:hbottom]+kwargs[:hbeam]/2],preserveframe=true)
        if i==1
            #scrootch a little bit back to make the overlap right
            seg=translate(seg,[-kwargs[:overlap]/4,0u"µm",0u"µm"],preserveframe=true)
        end
        return seg
    end
    #the center of the beam in xy plane
    cbeam = [(startx+stopx)/2,0u"µm"]
    #we can make the right side of the bridge by rotating leftsegs around cbeam
    rightsegs = map(leftsegs) do ls
        rotate(ls,pi,cbeam,preserveframe=true)
    end
    #the keystone is cut differently
    keystone = box(lkey,width,kwargs[:hbeam],kwargs[:dslice],
                   chamfer =[-kwargs[:cutangle]   -kwargs[:cutangle]
                             kwargs[:cutangle]    kwargs[:cutangle]])
    #move it into position
    keystone = translate(keystone,
                         vcat(cbeam,kwargs[:hbottom]+kwargs[:hbeam]/2),
                         preserveframe=true)
    #return a namedtuple so we can be fancy about how we write these
    Beam(leftsegs,rightsegs,keystone)
end

function Tessen.translate(b::Beam,args...;kwargs...)
    Beam([translate(x,args...;kwargs...) for x in b.leftsegs],
         [translate(x,args...;kwargs...) for x in b.rightsegs],
         translate(b.keystone,args...;kwargs...))
end

function Tessen.rotate(b::Beam,args...;kwargs...)
    Beam([rotate(x,args...;kwargs...) for x in b.leftsegs],
         [rotate(x,args...;kwargs...) for x in b.rightsegs],
         rotate(b.keystone,args...;kwargs...))
end


"""
```julia
post(;kwargs...)
```
Build a `Block` representing a square post
"""
function post(;kwargs...)
    #the z coordinates of our slices
    zcoords = range(start=zero(kwargs[:hbottom]),
                    stop=kwargs[:hbottom]+kwargs[:hbeam],
                    step=kwargs[:dslice])
    #`lengthpairs` will be a z => sidelength iterator
    lengthpairs = map(zcoords) do z
        z => kwargs[:wpost] - (if z < kwargs[:hbottom]
                                   #we're underneath the beams
                                   #the amount of undercut
                                   (kwargs[:hbottom]-z)*tan(kwargs[:chamferbottom])
                               else
                                   #we're on top of the beams
                                   #the amount of 'overcut'
                                   (z-kwargs[:hbottom])*tan(kwargs[:chamfertop])
                               end)
    end

    slicepairs = map(lengthpairs) do (z,l)
        coords = l/2 .* [[1,1],[-1,1],[-1,-1],[1,-1]]
        z => Slice([polycontour(coords)])
    end
    Block(slicepairs...)
end

struct Kernel
    support
    hammock
end

"""
```julia
kernel(;px,py,knx,kny,left=false,right=false,top=false,bottom=false,kwargs...)
```
Build one kernel, the directional arguments specify which sides get dangling
beams. The origin of the resulting `SuperBlock` will be at the centroid of the posts
"""
function kernel(;px,py,knx,kny,left=false,right=false,top=false,bottom=false,kwargs...)
    #write the top left post at (0,0), we'll scrooch everything around later
    firstpost = post(;kwargs...)
    posts = map([(i,j) for i in 0:(knx-1) for j in 0:(kny-1)]) do (i,j)
        #scrooch into position
        translate(firstpost,[i*px,j*py],preserveframe=true)
    end
    #merge the posts into one block
    postblock = merge(posts...)
    #the first horizontal beam is centered on (px/2,0) and has length `px-wpost`
    lhbeam = px - kwargs[:wpost]
    nhsegs = ceil(Int,lhbeam / kwargs[:maxseglength])
    if iseven(nhsegs)
        nhsegs += 1
    end
    #the first vertical beam is centered on (0,py/2) and has length `py-wpost`
    lvbeam = py - kwargs[:wpost]
    nvsegs = ceil(Int,lvbeam / kwargs[:maxseglength])
    if iseven(nvsegs)
        nvsegs += 1
    end
    #we'll give all of the beams the same number of segments for simplicity
    nsegs = max(nhsegs,nvsegs)
    firsthbeam = Beam(nsegs,kwargs[:wpost]/2,px-kwargs[:wpost]/2,kwargs[:wpost];kwargs...)
    #make one row's worth of the horizontal beams
    hbeamrow = map(0:(knx-2)) do i
        translate(firsthbeam,[i*px,zero(px)],preserveframe=true)
    end
    #do we have left beams
    if left
        pushfirst!(hbeamrow,translate(firsthbeam,[-px,zero(px)],preserveframe=true))
    end
    #do we have right beams
    if right
        push!(hbeamrow,translate(hbeamrow[end],[px,zero(px)],preserveframe=true))
    end
    #copy this row to make the others
    hbeams = vcat(([translate(hr,[zero(py),py*j],preserveframe=true)
                    for hr in hbeamrow] for j in 0:(kny-1))...)

    #same approach for the vertical beams
    firstvbeam = Beam(nsegs,kwargs[:wpost]/2,py-kwargs[:wpost]/2,kwargs[:wpost];kwargs...)
    #this beam is currently horizontal, do a rotation
    firstvbeam = rotate(firstvbeam,pi/2,preserveframe=true)
    #make one column's worth of vertical beams
    vbeamcol = map(0:(kny-2)) do j
        translate(firstvbeam,[zero(py),j*py],preserveframe=true)
    end
    #do we have top beams
    if top
        pushfirst!(vbeamcol,translate(firstvbeam,[zero(py),-py],preserveframe=true))
    end
    #do we have bottom beams
    if bottom
        push!(vbeamcol,translate(vbeamcol[end],[zero(py),py],preserveframe=true))
    end
    #copy this col to make the others
    vbeams = vcat(([translate(vb,[px*i,zero(px)],preserveframe=true)
                    for vb in vbeamcol] for i in 0:(knx-1))...)
    #we now want to merge all of our beams so the segments print starting from
    #the posts
    #get a vector containing each 'half' beam
    halfbeams = vcat(map([hbeams,vbeams]) do somebeams
                         #this vcat will result in a vector of vector where each entry is
                         #a half beam from `somebeams`
                         vcat(map(somebeams) do b
                                  [b.leftsegs, b.rightsegs]
                              end...)
                     end...)
    halfbeamblocks = map(zip(halfbeams...)) do hb
        #hb should now be all of the first segs, or all of the second segs, etc
        merge(hb...)
    end
    #now need to collect all of the keystones
    allkeys = vcat(map([hbeams,vbeams]) do somebeams
                       map(somebeams) do b
                           b.keystone
                       end
                   end...)
    #and merge them
    keyblock = merge(allkeys...)

    #now build the hammocks
    #figure out the total height of all of our layers
    hham = kwargs[:dhammockslice]*(kwargs[:nhammock]-1)
    @assert hham < kwargs[:hbeam] "hammocks can not be thicker than beams"
    #inital `overcut` due to the chamfered beams
    ocut_init = (kwargs[:hbeam] - hham) * tan(kwargs[:chamferbottom])
    #additional overcut per slice
    ocut_slice = kwargs[:dhammockslice] * tan(kwargs[:chamferbottom])
    #vertices of the first hammock if it was at the bottom of the beams, with no overlap
    bottom_verts = [
        [kwargs[:wpost]/2,-kwargs[:wpost]/2],
        [px-kwargs[:wpost]/2,-kwargs[:wpost]/2],
        [px-kwargs[:wpost]/2,py-kwargs[:wpost]/2],
        [kwargs[:wpost]/2,py-kwargs[:wpost]/2]
    ]
    #what sign should the cuts have at each vertex to dialate the hammock correctly
    dia_sign = [
        [-1,1], #top left vertex should move left and up
        [1,1], #top right vertex should move right and up
        [1,-1], #bottom right vertex should move right and down
        [-1,-1] #bottom left vertex should move left and down
    ]
    #first set ov vert coords
    verts = bottom_verts + (ocut_init + kwargs[:overlap])*dia_sign
    #initial z
    z = kwargs[:hbottom] + (kwargs[:hbeam] - hham)
    hamslices = []
    for i in 1:kwargs[:nhammock]
        #add this slice to our vector
        push!(hamslices, z=>Slice([polycontour(verts)]))
        #update for the next round
        z += kwargs[:dhammockslice]
        verts += ocut_slice*dia_sign
    end
    #this is a hammock to the bottom right of the first post
    firsthamblock = Block(hamslices...)
    #make the first row of hammocks
    hamrow = map(0:(knx - 2)) do i
        translate(firsthamblock,[i*px,zero(px)],preserveframe=true)
    end
    #add hammocks for dangling beams
    if left
        pushfirst!(hamrow,translate(firsthamblock,[-px,zero(px)],preserveframe=true))
    end
    if right
        push!(hamrow,translate(hamrow[end],[px,zero(px)],preserveframe=true))
    end

    #now copy this row
    hams = map(0:(kny-2)) do j
        [translate(hr,[zero(py),j*py],preserveframe=true) for hr in hamrow]
    end

    #add on if we have top or bottom beams
    if top
        pushfirst!(hams,[translate(hr,[zero(py),-py],preserveframe=true) for hr in hamrow])
    end
    if bottom
        push!(hams,[translate(he,[zero(py),py],preserveframe=true) for he in hams[end]])
    end

    #merge all the hammocks into one block
    hamblock = merge(vcat(hams...)...)
    #hatch the hamblock
    @info "hatching a kernel"
    hatchedham = hatch(hamblock;dhatch=kwargs[:dhammockhatch],
                       bottomdir = pi/4,
                       diroffset = pi/2)

    #build a support `SuperBlock`
    sb = SuperBlock(postblock,halfbeamblocks...,keyblock)
    
    #move the origin to the center of the posts
    neworigin = [(knx-1)*px/2,(kny-1)*py/2]
    sbt = translateorigin(sb,neworigin)
    ht = translateorigin(hatchedham,neworigin)
    #hatch and return
    sbth = hatch(sbt,dhatch=kwargs[:dhatch],
                 bottomdir = pi/4,
                 diroffset = pi/2)
    Kernel(sbth,ht)
end

"""
```julia
scaffold(scaffolddir[, configfilename])
scaffold(scaffolddir,configdict)
```
Build a directory of .gwl files to build a scaffold in a directory at `scaffolddir`. The
scaffold's geometrical parameters can be provided as a `Dict` or read from `configfilename`.
If `configfilename` is omitted, the parameters will be read from `"config.jl"`
"""
function scaffold end

function scaffold(scaffolddir,kwargs::Dict)
    #make a folder to hold the files for this scaffold
    mkdir(scaffolddir)
    #move down in there do compile the geometry
    compgeom = cd(scaffolddir) do
        #make a folder for our pre-compiled scaffold parts
        mkdir("scripts")
        #assemble and pre-compile the bumpers
        b=bumper(;kwargs...)
        #this bumper is currently at the origin, translate it up into position
        btop=translate(b,[0u"µm",(kwargs[:wscaf]-kwargs[:wbumper])/2,0u"µm"])
        #get the bottom by rotation
        bbottom = rotate(btop,pi)
        @info "compiling bumpers"
        tophatched=hatch(btop,dhatch=kwargs[:dhatch],bottomdir=pi/4)
        bottomhatched=hatch(bbottom,dhatch=kwargs[:dhatch],bottomdir=pi/4)
        #printing the 'bottom' bumper first reduces stage movement
        bumpers = CompiledGeometry(joinpath("scripts","bumpers.gwl"),bottomhatched,tophatched;
                                   laserpower=kwargs[:laserpower],
                                   scanspeed=kwargs[:scanspeed])
        #determine our kernel parameters
        #maximum possible pitch
        pmax = kwargs[:wpost] + kwargs[:lbeammax]
        #calculate our grid dimensions
        
        ny = ceil(Int,
                  #distance from the top of the top row of posts to the bottom bumper
                  (kwargs[:wscaf] - 2*kwargs[:wbumper] - kwargs[:lbeammax])/
                      pmax)
                nx = #start with one post
            1 + 
            ceil(Int,
                 #distance from the end of the leftmost post to the end of the scaffold
                 #we have to subtract off wbumper to account for the rounded ends of the
                 #bumper.
                 (kwargs[:lscaf] - kwargs[:wpost] - kwargs[:wbumper])/
                     pmax)

        #we don't want nx or ny to be prime
        (nx,ny) = map([nx,ny]) do ni
            isprime(ni) ? ni+1 : ni
        end

        #get the actual beam lengths
        nbeamx = nx - 1
        nbeamy = ny + 1
        
        lbeamx = (kwargs[:lscaf] - nx*kwargs[:wpost] - kwargs[:wbumper])/nbeamx
        lbeamy = (kwargs[:wscaf] - 2*kwargs[:wbumper] - ny*kwargs[:wpost])/nbeamy
        @assert lbeamx <= kwargs[:lbeammax]
        @assert lbeamy <= kwargs[:lbeammax]

        #calculate actual pitch
        px = lbeamx + kwargs[:wpost]
        py = lbeamy + kwargs[:wpost]
        #figure out the biggest kernel which fits neatly into our grid
        #define a helper function so we can get all multiples of our nx and ny
        function allmultiples(x)
	    filter(1:x) do y
	        (x%y) == 0
	    end
        end

        #what is the biggest 'kernel' of posts that 1) fits evenly into the desired space
        #and can fit inside a circle with diameter `dfield`
        #we will assume each 'kernel' includes connecting bridges on one side, the top and bottom.
        #the bottom bridges will only actually be written on the bottom row

        allnumcombos = [(nxi,nyi) for nxi in allmultiples(nx), nyi in allmultiples(ny)]
        possiblenumcombos = filter(allnumcombos) do (nxi,nyi)
	    sizex = nxi*px + lbeamx + 2*kwargs[:overlap] #two 'dangling' beams
	    sizey = nyi*py + lbeamy + 2*kwargs[:overlap] #two 'dangling' beams
	    sqrt(sizex^2 + sizey^2) < kwargs[:dfield] #does this fit
        end
        #get the option which contains the most posts
        (_,maxcomboindex) = map(possiblenumcombos) do (nx,ny)
	    nx*ny
        end |> findmax
        #assign this max size to variables for convenience
        (knx,kny) = possiblenumcombos[maxcomboindex]
        
        #we will write kernels in a serpentine pattern, top to bottom, connecting
        #to kernels 'above' and 'behind' us with beams. On the last row we will connect
        #to the bumper below us

        #get the coordinates of the center of each kernel.
        #first get the size of the whole shebang
        postarraydims = [knx*px + kwargs[:wpost], kny*py + kwargs[:wpost]]
        #the array of posts is centered on (0,0)
        arraycorner = [-postarraydims[1]/2,postarraydims[1]/2]
        #coordinates of the first kernel center
        firstkernelcenter = arraycorner + [
            kwargs[:wpost]/2 + (knx-1)*px/2,
            -kwargs[:wpost]/2 - (kny-1)*py/2
        ]
        #compile all the kernels we are going to need
        common_args = Dict(:px => px, :py => py, :knx => knx, :kny => kny)
        kernelflavors = [
            #kernel(;px,py,knx,kny,left=false,right=false,top=false,bottom=false,kwargs...)
            :t => kernel(;top=true, common_args...,kwargs...),
            :lt => kernel(;top=true,left=true, common_args...,kwargs...),
            :rt => kernel(;top=true,right=true, common_args...,kwargs...),
            :tb => kernel(;top=true,bottom=true, common_args...,kwargs...),
            :ltb => kernel(;left=true, top=true,bottom=true, common_args...,kwargs...),
            :rtb => kernel(;right=true, top=true,bottom=true, common_args...,kwargs...)
        ]
        #compile all our flavors and store in a dict
        kernels = map(kernelflavors) do (k,v)
            #need to change how our kernel function works so we can change laser power
            #for the hammocks
            @info "Compiling $v kernel"
            k => [CompiledGeometry("$(k)_support.gwl",v.support,laserpower=kwargs[:laserpower],
                                   scanspeed=kwargs[:scanspeed]),
                  CompiledGeometry("$(k)_hammock.gwl",v.hammock,laserpower=kwargs[:hamlaserpower],
                                   scanspeed=kwargs[:hamscanspeed])
             ]            
        end |> Dict
        #need to vcat since each kernel has support and a hammock
        allkernels = []
        i = 0
        j = 0
        while true
            while true
                kernelcenter = firstkernelcenter + [knx*px*i,-kny*py*j]
                if ((i==0) && iseven(j)) || ((i == (knx-1)) && isodd(j))
                    #first kernel on a new row, just want beams on top,
                    #unless this is the last row
                    thiskernel = j == kny-1 ? kernels[:tb] : kernels[:t]
                    thiskerneltranslated = [translate(k,kernelcenter) for k in thiskernel]
                    push!(allkernels,thiskerneltranslated...)
                else
                    #not the first kernel, want beams on top and on the side
                    thiskernel = if iseven(j)
                        #beams on the left
                        j == kny-1 ? kernels[:ltb] : kernels[:lt]
                    else
                        #beams on the right
                        j == kny-1 ? kernels[:rtb] : kernels[:rt]
                    end
                    thiskerneltranslated = [translate(k,kernelcenter) for k in thiskernel]
                    push!(allkernels,thiskerneltranslated...)
                end
                if iseven(j)
                    i+=1
                else
                    i-=1
                end
                if i>(knx-1)
                    break
                end
            end
            j+=1
            if j>(kny-1)
                break
            end
        end
        vcat(bumpers,allkernels)
    end
    GWLJob(joinpath(scaffolddir,"scaffold.gwl"),compgeom...;
           stagespeed=kwargs[:stagespeed],interfacepos=kwargs[:interfacepos])
end

function scaffold(scaffolddir,configfilename::String)
    config = include(configfilename)
    scaffold(scaffolddir,config)
end

scaffold(scaffolddir) = scaffold(scaffolddir,"config.jl")

end # module Pancake
