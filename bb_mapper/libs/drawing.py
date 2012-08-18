from reportlab.pdfgen import canvas
from reportlab.lib.units import cm
from reportlab.lib.colors import black, white, HexColor
from loaders import load_genbank
from config import fct_flags, fct_colors, idpt, min_size
from array_tetris import offset_coord, nudge_coord, shade_split, \
    coord_flipper

### Hardcoded presets ###
## General proportions ##
u 		=  0.01 #1/50 	# conversion factor for sequence length-related values
## (0.0125, increase for short sequences)
hmar 	= 2*cm			# horizontal margin to canvas (1*cm)
vmar 	= 2.5*cm		# vertical margin to canvas (2*cm)
pNsize 	= 5*cm			# width to set aside for plasmid names (2.5*cm)
di 		= 0.18*cm		# half-length of interruption/frameshift tick marks (0.18*cm)
doLup 	= 0.6*cm		# distance of ORF labels from their ORF (0.6*cm) (above)
doLdn   = 0.8*cm		# distance of ORF labels from their ORF (0.8*cm) (below)
minL 	= 1000			# minimum ORF size (in base pairs) for 'full arrows'
                        # (700, decrease for short seq)
w 		= 0.3*cm     	# half-width of the tail (0.3*cm)
h 		= 0.175*cm		# distance from the side tips of the head to the neck (0.175*cm)
dBL	 	= 3*cm			# distance between plasmid baselines (4*cm)
da 		= 0.4*cm		# distance of alignment tick marks from the corresponding plasmid baseline (0.4*cm)
ck_hsp 	= 170		    # horizontal spacing between columns (170)
ck_vsp 	= 17		    # vertical spacing between items (17)
tm 		= 0.2*cm		# size of alignment tick marks (0.2*cm)
y_adj   = -4.2*cm       # vertical adjustment factor

## Typeface, fonts, sizes ##
rFont 	= "Helvetica"			# regular type ("Helvetica")
bFont 	= "Helvetica-Bold"		# bold type ("Helvetica-Bold")
LfSize 	= 14					# large font size (14)
NfSize 	= 12					# normal font size (12)
SfSize 	= 10					# small font size (10)

## Special feature settings ##
SFX 	= 'off'			# allow drawing special sequence features beside CDS/ORF ('on' or 'off')
osym 	= '*'			# symbol for origin of sequence (*)
snp		= '*'			# symbol for SNP (*)

## Scale of sequence, in legend ##
scX 	= 0*u			# starting point on the X axis(0)
incrT 	= 1000*u		# increment size in basepairs (1000, signifying 1 kb)
incrN 	= 5				# increment number (5, or 1 for short sequences)
dip 	= -2*cm			# starting point on Y axis (-2*cm)
dop 	= 0.3*cm		# size of vertical tick

## Color key settings ##
ckX 		= 5*cm		# starting point on X axis (5*cm)
ckY 		= dip		# starting point on Y axis (dip, same as legend scale, but can take any numerical value)
ck_hsp 		= 170		# horizontal spacing between columns (170)
ck_vsp 		= 17		# vertical spacing between items (17)
ck_boxX 	= 12		# horizontal side size of box (10)
ck_boxY 	= 12		# vertical side size of box (10)
ck_htxoff 	= 20		# horizontal offset for text from box origin (1)
ck_vtxoff 	= 1			# vertical offset for text from box origin (20)
lay_MAX 	= 2			# maximum items in a column (4, can also be a fraction of len(functions))

## Heatmap key settings ##
hk_u		= 3		    # length conversion factor (1.5)
hkY			= dip		# starting point on Y axis (dip)
hk_boxX		= ck_boxX*2	# horizontal side size of box (ck_boxY)


def canvasser(hCan,vCan,transX,transY,outfile) :
    """Initialize canvas."""
    canvasN = canvas.Canvas(outfile, pagesize=(hCan,vCan))
    canvasN.translate(transX,transY)
    canvasN.setStrokeColor(black)
    canvasN.setFillColor(white)
    canvasN.setLineWidth(1)
    canvasN.setLineJoin(1)
    canvasN.setLineCap(0)
    return canvasN

def base_draw(canvas, genome, feats, key, dop_Y, Y0, X_shift, map_mode,
             annot_cnt, seq_len, annot_mode, side) :
    """Draw contig baseline and features."""
    # unpack info
    cName = genome.name
    cLen = genome.len
    nudge = genome.nudge
    offset = genome.offset
    # draw plasmid baseline
    baseliner(cLen, canvas, Y0, nudge)
    # label the baseline with plasmid name and size
    labeller(cName, cLen, canvas, Y0)
    # label the annotations list
    Y_annot = 0
    if map_mode != 'n':
        canvas.setFont(bFont, LfSize)
        canvas.drawString(X_shift, y_adj, cName)
    canvas.setFont(rFont, SfSize)
    Y_annot +=2
    # filter and draw annotation features
    ORFcnt = 0
    shift_flag = False
    for feature in feats :
        if feature.type == 'contig':
            contig_ticker(canvas, feature, cLen, Y0, nudge, offset, side)
        elif feature.type == 'spacer':
            spacer_ticker(canvas, feature, cLen, Y0, nudge, offset, side)
        elif feature.type == 'ref_seg':
            ref_ticker(canvas, feature, cLen, Y0, nudge, offset)
        elif feature.type == 'CDS' or feature.type == 'cds':
            ORFcnt += 1
            # determine functional category color
            try:
                annot = feature.qualifiers.get(key)[0]
            except TypeError:
                annot = 'none'
            fct_key = annot_color(fct_flags, annot)
            color_hex = HexColor(fct_colors[fct_key][0])
            # calculate coordinates for canvas
            featL, midLZ, coords, split_flag = orf_coords(feature, Y0, cLen,
                                                          nudge, offset)
            # check for CDS sitting across the origin
            if split_flag:
                # split coords
                coords1, featL1, coords2, featL2 = orf_split(coords, cLen)
                # draw square feature
                orf_eus(canvas, featL1, coords1, color_hex, shape='square')
                # draw arrow feature
                orf_eus(canvas, featL2, coords2, color_hex, shape=None)
            else:
                # draw arrow feature
                orf_eus(canvas, featL, coords, color_hex, shape=None)
            # write annotation line and CDS number
            if map_mode != 'n':
                if map_mode == 'single' and Y_annot-1 > annot_cnt/2:
                    X_shift = seq_len/2
                    if not shift_flag:
                        Y_annot -= annot_cnt/2
                        shift_flag = True
                Y_annot = orf_annot(canvas, cLen, annot, Y_annot, ORFcnt,
                                   Y0+dop_Y, midLZ, X_shift, split_flag,
                                   annot_mode)

def baseliner(cLen, canvas, Y_map, nudge):  # loop offset has no effect
    """Draw sequence baseline."""
    canvas.setLineWidth(3)
    y0 = Y_map
    Zs = 0              # all sequences are initially aligned on the left
    Ze = cLen
    # calculate nudge offsets
    offZs = nudge_coord(Zs, nudge)
    offZe = nudge_coord(Ze, nudge)
    # set zeroes
    x0 = offZs*u
    x1 = offZe*u
    pBL = canvas.beginPath()
    pBL.moveTo(x0,y0)
    pBL.lineTo(x1,y0)
    canvas.drawPath(pBL, stroke=1, fill=0)
    pBL.close()

def labeller(cName, cLen, canvas, Y_map) :
    """Label baselines with genome/contig name and size."""
    canvas.setFillColor(black)
    y0 = Y_map
    x0 = -pNsize                # retreat into the left margin to write out name and size
    y1 = y0 + ck_vsp/10         # bump name up a bit from BL level
    y2 = y0 - ck_vsp            # bump size down a bit from BL level
    pLenStr = str(float(cLen)/1000) # converts number to suitable form
    canvas.setFont(bFont,LfSize)
    canvas.drawString(x0,y1,cName)
    canvas.setFont(rFont,NfSize)
    canvas.drawString(x0,y2,pLenStr+' kb')

def orf_annot(canvas, cLen, annot, Y_annot, ORFcnt, cnt_Y, midLZ, X_shift,
             split_flag, annot_mode):
    """Write annotation to feature list."""
    flag = False
    if not annot_mode == 'all':
        if annot == 'no match' or annot == 'hypothetical protein' :
            flag = True
    if not flag:
        # write CDS numbers
        if split_flag:
            canvas.drawCentredString(1*u, cnt_Y, "]"+str(ORFcnt))
            canvas.drawCentredString(cLen*u, cnt_Y, str(ORFcnt)+"[")
        else:
            canvas.drawCentredString(midLZ, cnt_Y, str(ORFcnt))
        # write annotation
        y_annot_adj = y_adj-(Y_annot*ck_vsp)
        canvas.setFont(rFont, SfSize)
        canvas.drawString(X_shift, y_annot_adj, str(ORFcnt)+'. '+annot)
        canvas.setFont(rFont, NfSize)
        Y_annot +=1
    return Y_annot

def orf_split(coords, cLen):
    """Split CDS that sit across the map origin."""
    xs, xe, xn, y0, yt, yb, ynt, ynb = coords
    if xs < xe:
        coords1 = xs, 1*u, xn, y0, yt, yb, ynt, ynb
        coords2 = cLen*u, xe, xn, y0, yt, yb, ynt, ynb
        featL1 = xs/u
        featL2 = cLen-xe/u
    else:
        coords1 = xs, cLen*u, xn, y0, yt, yb, ynt, ynb
        coords2 = 1*u, xe, xn, y0, yt, yb, ynt, ynb
        featL1 = cLen-xs/u
        featL2 = xe/u
    return coords1, featL1, coords2, featL2

def orf_coords(feature, Y_map, cLen, nudge, offset):
    """Calculate CDS coordinates in drawing space."""
    # evaluate what strand the ORF is on (determines direction of arrow)
    cstrand = feature.strand
    if cstrand is None:
        cstrand = 1
    # take start and end points
    location = feature.location
    Zs = location.nofuzzy_start
    Ze = location.nofuzzy_end
    featL = Ze - Zs
    # calculate loop offset coordinates
    loop_offZs = offset_coord(Zs, cLen, offset)
    loop_offZe = offset_coord(Ze, cLen, offset)
    # calculate nudge offset coordinates
    offZs = nudge_coord(loop_offZs, nudge)
    offZe = nudge_coord(loop_offZe, nudge)
    # calculate X axis coordinates (expr of cstrand has changed)
    if cstrand == -1 :	# reverse orientation
        xs,xe = offZe*u,offZs*u		# start and end
        xn = xe+minL*u		        # neck of arrow
    else :				# forward orientation
        xs,xe = offZs*u,offZe*u		# start and end
        xn = xe-minL*u		        # neck of arrow
    midLZ = ((offZs+offZe)/2)*u	    # middle of ORF for optional label
    # evaluate splits
    split_flag = False
    if (xs < xe and cstrand == -1) or (xs > xe and cstrand == 1):
        midLZ = 0
        split_flag = True
    # set Y axis coordinates
    y0 = Y_map
    yt,yb,ynt,ynb = y0+w,y0-w,y0+h,y0-h
    coords = xs, xe, xn, y0, yt, yb, ynt, ynb
    return featL, midLZ, coords, split_flag

def orf_eus(canvas, featL, coords, color_hex, shape):
    """Draw CDS and write count."""
    xs, xe, xn, y0, yt, yb, ynt, ynb = coords
    canvas.setLineWidth(1)
    # initialize path
    pORF = canvas.beginPath()
    if shape == 'square':
        pORF.moveTo(xs,ynt)
        pORF.lineTo(xe,ynt)
        pORF.lineTo(xe,ynb)
        pORF.lineTo(xs,ynb)
        pORF.lineTo(xs,ynt)
    # draw triangle-shaped ORFS
    elif featL <= minL:
        pORF.moveTo(xs,yt)
        pORF.lineTo(xe,y0)
        pORF.lineTo(xs,yb)
        pORF.lineTo(xs,yt)
    # draw arrow-shaped ORFS
    else:
        pORF.moveTo(xs,ynt)
        pORF.lineTo(xn,ynt)
        pORF.lineTo(xn,yt)
        pORF.lineTo(xe,y0)
        pORF.lineTo(xn,yb)
        pORF.lineTo(xn,ynb)
        pORF.lineTo(xs,ynb)
        pORF.lineTo(xs,ynt)
    # evaluate function category and set fill color
    canvas.setFillColor(color_hex)
    # finalize object path
    canvas.drawPath(pORF, stroke=1, fill=1)
    pORF.close()
    canvas.setFillColor(black)

def contig_ticker(canvas, feature, cLen, Y0, nudge, offset, side):
    """Draw contig separators."""
    # get contig name
    name = feature.qualifiers.get('locus_tag')[0]
    # take start and end points
    location = feature.location
    Zs = location.nofuzzy_start
    Ze = location.nofuzzy_end
    # calculate loop offset coordinates
    loop_offZs = offset_coord(Zs, cLen, offset)
    loop_offZe = offset_coord(Ze, cLen, offset)
    # calculate nudge offset coordinates
    offZs = nudge_coord(loop_offZs, nudge)
    offZe = nudge_coord(loop_offZe, nudge)
    xs, xe = offZs*u, offZe*u
    # set Y axis coordinates
    if side == 'low':
        y0 = Y0-dop*3.5
        yT1 = y0-h*4
        yT2 = y0-h*4.5
        yT3 = y0-h*8
    else:
        y0 = Y0+dop*3.5
        yT1 = y0+h*4
        yT2 = y0+h*4.5
        yT3 = y0+h*8
    # draw
    canvas.setLineWidth(2)
    ttl = canvas.beginPath()
    ttl.moveTo(xs,y0)
    ttl.lineTo(xs,yT1)
    ttl.lineTo(xs+dop,yT1)
    canvas.drawPath(ttl, stroke=1, fill=0)
    canvas.setFont(bFont, NfSize)
    canvas.drawString(xs+dop*2, yT2, name)
    canvas.setFont(rFont, SfSize)
    canvas.drawString(xs+dop*2,yT3,"".join(["[",str(Zs),"-",str(Ze),"]"]))
    canvas.setFont(rFont, NfSize)
    ttl.close()

def spacer_ticker(canvas, feature, cLen, Y0, nudge, offset, side):
    """Draw separator indicators."""
    # take start and end points
    location = feature.location
    Zs = location.nofuzzy_start
    Ze = location.nofuzzy_end
    # calculate loop offset coordinates
    loop_offZs = offset_coord(Zs, cLen, offset)
    loop_offZe = offset_coord(Ze, cLen, offset)
    # calculate nudge offset coordinates
    offZs = nudge_coord(loop_offZs, nudge)
    offZe = nudge_coord(loop_offZe, nudge)
    xs, xe = offZs*u, offZe*u
    # set Y axis coordinates
    y0 = Y0-dop*3.5
    if side == 'low':
        yT1 = y0-h*4
    else:
        yT1 = y0+h*4
    # draw
    canvas.setLineWidth(2)
    ttl = canvas.beginPath()
    ttl.moveTo(xs,y0)
    ttl.lineTo(xs,yT1)
    ttl.lineTo(xe,yT1)
    ttl.lineTo(xe,y0)
    canvas.drawPath(ttl, stroke=1, fill=0)
    ttl.close()

def ref_ticker(canvas, feature, cLen, Y0, nudge, offset):
    """Draw contig separators."""
    # get contig name
    name = feature.qualifiers.get('id')[0]
    # take start and end points
    location = feature.location
    Zs = location.nofuzzy_start
    Ze = location.nofuzzy_end
    # calculate loop offset coordinates
    loop_offZs = offset_coord(Zs, cLen, offset)
    loop_offZe = offset_coord(Ze, cLen, offset)
    # calculate nudge offset coordinates
    offZs = nudge_coord(loop_offZs, nudge)
    offZe = nudge_coord(loop_offZe, nudge)
    xs, xe = offZs*u, offZe*u
    xmid = (xe+xs)/2
    # set Y axis coordinates
    y0 = Y0+dop*3
    # draw
    canvas.setLineWidth(2)
    ttl = canvas.beginPath()
    ttl.moveTo(xs,y0+w)
    ttl.lineTo(xs,y0+w+h*2)
    ttl.lineTo(xe,y0+w+h*2)
    ttl.lineTo(xe,y0+w)
    canvas.drawPath(ttl, stroke=1, fill=0)
    canvas.setFont(bFont, NfSize)
    canvas.drawCentredString(xmid, y0+h*5, name)
    canvas.setFont(rFont, NfSize)
    ttl.close()

def seq_scale(canvas, scX, scY, incrT, incrN, dip, dop) :
    """Draws the sequence scale bar."""
    canvas.setLineWidth(1.2)
    canvas.setFillColor(black)
    incrCNT = 0							    # initialize count of increments
    psc = canvas.beginPath()
    psc.moveTo(scX, scY+dip-dop)		        # start at beginning (duh!)
    psc.lineTo(scX+incrT*incrN, scY+dip-dop)	# draw the scale baseline
    # draw ticks until the max number of increments is reached
    while incrCNT <= incrN :
        psc.moveTo(scX+incrT*incrCNT, scY+dip-dop)
        psc.lineTo(scX+incrT*incrCNT, scY+dip)
        incrCNT += 1
    canvas.drawPath(psc, stroke=1, fill=0)
    psc.close()
    # write out scale extremities values (needs hand-fix if not using kbs)
    canvas.setFont(rFont,NfSize)
    canvas.drawRightString(scX, scY+dip+dop, '0')
    canvas.drawString(scX+incrT*incrN, scY+dip+dop, str(incrN)+' kb')

def annot_color(fct_flags, annotation):
    """Look up the color to use based on annotation keywords."""
    annot_line = annotation.lower()
    fct_key = 'oth'
    for key in fct_flags:
        i = 0
        while i < len(fct_flags[key]):
            if annot_line.find(fct_flags[key][i]) > -1:
                fct_key = key
                break
            i +=1
    return fct_key

def contig_draw(contig, in_file, out_file, annot_mode, key):
    """Draw sequence map of a single contig to file."""
    # load contig record
    seq_record = load_genbank(in_file)
    ctg_len = len(seq_record.seq)
    feats = seq_record.features
    cds = [feature for feature in feats
           if feature.type == 'CDS' or feature.type == 'cds']
    if annot_mode == 'all':
        annot_cds = [len(cds)]
    else:
        try:
            annot_cds = [1 for feature in cds
                         if feature.qualifiers.get(key)[0] != 'no match']
        except TypeError:
            annot_cds = []
    annot_cnt = sum(annot_cds)
    # calculate main canvas dimensions
    if ctg_len*u < 2000:
        seq_len = 2000
    else:
        seq_len = ctg_len*u
    hCan = hmar*2 + pNsize + seq_len
    vCan = dBL + vmar*4 + (annot_cnt/2)*ck_vsp
    transX = hmar + pNsize
    transY = dBL + vmar*2 + (annot_cnt/2)*ck_vsp
    ctg_Y = vmar
    # set up main canvas
    canvas = canvasser(hCan, vCan, transX, transY, out_file)
    # draw contig baseline and features -- offsetting disabled
    base_draw(canvas, contig, feats, key, -doLdn, ctg_Y, 0, 'single',
             annot_cnt, seq_len, annot_mode, 'top')
    # draw scale
    seq_scale(canvas, (ctg_len*u)-pNsize, incrT, incrN, dip, dop)
    # write to file and finalize the figure
    canvas.showPage()
    canvas.save()

def pairwise_draw(ref, query, segs, map_file, mode1, mode2, annot_mode,
                  key1, key2):
    """Draw pairwise alignment map with similarity shading."""
    # load ref record
    ref_record = load_genbank(ref.gbk)
    ref_feat = ref_record.features
    ref_cds = [feature for feature in ref_feat
               if feature.type == 'CDS' or feature.type == 'cds']
    if annot_mode != 'all':
        try:
            ref_annot_cds = [1 for cds in ref_cds
                             if cds.qualifiers.get(key1)[0] !=
                                'hypothetical protein' and \
                                cds.qualifiers.get(key1)[0] !=
                                'no match']
        except TypeError:
            ref_annot_cds = []
        ref_annot_cnt = sum(ref_annot_cds)
    else:
        ref_annot_cnt = len(ref_cds)
    # load query record
    query_record = load_genbank(query.gbk)
    if query.invert:
        query_record = query_record.reverse_complement()
    q_feat = query_record.features
    query_cds = [feature for feature in q_feat
                 if feature.type == 'CDS' or feature.type == 'cds']
    if annot_mode != 'all':
        try:
            query_annot_cds = [1 for cds in query_cds
                               if cds.qualifiers.get(key2)[0] !=
                                'hypothetical protein' and \
                                cds.qualifiers.get(key2)[0] !=
                                'no match']
        except TypeError:
            query_annot_cds = []
        query_annot_cnt = sum(query_annot_cds)
    else:
        query_annot_cnt = len(query_cds)
    # calculate main canvas dimensions - horizontal
    if ref.len+ref.nudge > query.len:
        ctg_len = ref.len+ref.nudge
    else:
        ctg_len = query.len
    if ctg_len*u < 2000:
        seq_len = 2000
    else:
        seq_len = ctg_len*u
    hCan = hmar*2 + pNsize + seq_len
    # calculate main canvas dimensions - vertical
    if mode1 == 'single' and mode2 == 'n':
        annot_cnt = ref_annot_cnt
        annot_len = annot_cnt/2
    else:
        annot_cnt = max(ref_annot_cnt, query_annot_cnt)
        annot_len = annot_cnt
    vCan = dBL + vmar*6 + annot_len*ck_vsp
    transX = hmar + pNsize
    transY = dBL + vmar*1.8 + annot_len*ck_vsp
    ref_Y = vmar*2.8
    query_Y = vmar
    # set up main canvas
    m_canvas = canvasser(hCan, vCan, transX, transY, map_file)
    # draw scale
    seq_scale(m_canvas, (ctg_len*u)-pNsize, 0, incrT, incrN, dip, dop )
    # draw shading legend
    heatkey(m_canvas, -pNsize, -pNsize/2)
    # draw ref baseline and features
    base_draw(m_canvas, ref, ref_feat, key1, doLup, ref_Y, 0,
              mode1, annot_cnt, seq_len, annot_mode, 'top')
    # draw query baseline and features
    base_draw(m_canvas, query, q_feat, key2, -doLdn, query_Y, seq_len/2,
              mode2, annot_cnt, seq_len, annot_mode, 'low')
    # draw pairwise similarity shading
    try:
        for xa, xb, xc, xd, idp in segs:
            # evaluate color shading category
            sh_color = HexColor(simcolor(idp))
            # check for split
            if abs(xa) > abs(xb) or abs(xc) > abs(xd):
                new_segpairs = shade_split(xa, xb, xc, xd, ref, query)
                for xa1, xb1, xc1, xd1 in new_segpairs:
                    # draw shading
                    shadowfax(m_canvas, xa1, xb1, xc1, xd1, ref_Y, query_Y, sh_color)
            else:
                # draw shading
                shadowfax(m_canvas, xa, xb, xc, xd, ref_Y, query_Y, sh_color)
    except TypeError:
        pass
    # write to file and finalize the figure
    m_canvas.showPage()
    m_canvas.save()

def multi_draw(g_pairs, segdata_list, mapfile):
    """Draw multiple alignment map with similarity shading."""
    print "Lookin\' good!"
    # compile info
    lengths = [g_pairs[0][0].len+g_pairs[0][0].nudge]
    g_to_draw = [g_pairs[0][0]]
    for ref, query in g_pairs:
        lengths.append(query.len+query.nudge)
        g_to_draw.append(query)
    max_len = max(lengths)
    # calculate main canvas dimensions - horizontal
    if max_len*u < 2000:
        seq_len = 2000
    else:
        seq_len = max_len*u
    hCan = hmar*4 + pNsize + seq_len
    # calculate main canvas dimensions - vertical
    vCan = dBL*len(g_pairs) + vmar*4
    transX = hmar + pNsize
    transY = dBL*len(g_pairs)
    init_Y = vmar*2
    # set up main canvas
    m_canvas = canvasser(hCan, vCan, transX, transY, mapfile)
    # draw scale (max_len*u)-pNsize, hmar
    seq_scale(m_canvas, 2*hCan/3, -vmar*2, incrT, incrN, dip, dop)
    # draw shading legend
    heatkey(m_canvas, hCan-hmar*5, init_Y+vmar)
    # draw ref baseline and features
    counter = 0
    for genome in g_to_draw:
        g_record = load_genbank(genome.gbk)
        g_feat = g_record.features
        g_cds = [feature for feature in g_feat
                   if feature.type == 'CDS' or feature.type == 'cds']
        ref_Y = init_Y-dBL*counter
        base_draw(m_canvas, genome, g_cds, '', doLup, ref_Y, 0, 'n', 0,
                  seq_len, 'n', 'n')
        counter +=1
    counter = 0
    for ref, query in g_pairs:
        ref_Y = init_Y-dBL*counter
        query_Y = init_Y-dBL*(counter+1)
        # draw pairwise similarity shading
        try:                                  # TODO: adapt Y
            for xa, xb, xc, xd, idp in segdata_list[counter]:
                # evaluate color shading category
                sh_color = HexColor(simcolor(idp))
                # check for split
                if abs(xa) > abs(xb) or abs(xc) > abs(xd):
                    new_segpairs = shade_split(xa, xb, xc, xd, ref, query)
                    for xa1, xb1, xc1, xd1 in new_segpairs:
                        # draw shading
                        shadowfax(m_canvas, xa1, xb1, xc1, xd1, ref_Y,
                                  query_Y, sh_color)
                else:
                    # draw shading
                    shadowfax(m_canvas, xa, xb, xc, xd, ref_Y, query_Y,
                              sh_color)
            counter +=1
        except TypeError:
            pass
    # write to file and finalize the figure
    m_canvas.showPage()
    m_canvas.save()

def shadowfax(canvas_def, xa, xb, xc, xd, aby0, cdy0, sh_color):
    """Draw shaded area between homologous segments."""
    # cancel drawing if segments too small to draw
    if abs(xb)-abs(xa) < min_size:
        pass
    else:
        # convert sequence-scale values to canvas-space
        axr, bxr, cxr, dxr = xa*u, xb*u, xc*u, xd*u
        # check for negatives that need to be flipped
        ax, bx = coord_flipper(axr,bxr)
        cx, dx = coord_flipper(cxr,dxr)
        # these are the actual Y coordinates that will be used to draw the cues
        aby1 = aby0-da
        aby2 = aby1-tm
        cdy1 = cdy0+da
        cdy2 = cdy1+tm
        # this draws the parallelograms between matching segments
        canvas_def.setFillColor(sh_color)
        ppg = canvas_def.beginPath()
        ppg.moveTo(ax,aby2)
        ppg.lineTo(bx,aby2)
        ppg.lineTo(dx,cdy2)
        ppg.lineTo(cx,cdy2)
        ppg.lineTo(ax,aby2)
        canvas_def.drawPath(ppg, stroke=0, fill=1)
        ppg.close()
        # this draws the tick marks and lines delineating matching segments
        canvas_def.setLineWidth(0.2)
        puck = canvas_def.beginPath()
        puck.moveTo(ax,aby1)
        puck.lineTo(ax,aby2)
        puck.lineTo(cx,cdy2)
        puck.lineTo(cx,cdy1)
        puck.moveTo(bx,aby1)
        puck.lineTo(bx,aby2)
        puck.lineTo(dx,cdy2)
        puck.lineTo(dx,cdy1)
        canvas_def.drawPath(puck, stroke=1, fill=0)
        puck.close()

def simcolor(idp):
    """Evaluate class of similarity."""
    id_cats = [x for x in idpt if idp>=x]
    sh_hex = idpt[max(id_cats)]
    return sh_hex

def heatkey(canvas, hkX, hkY):
    """Draw color key for the heat map."""
    canvas.setLineWidth(1)
    canvas.setLineCap(0)
    canvas.setFillColor(black)
    canvas.setFont(bFont, LfSize)
    canvas.drawString(hkX, hkY, "Nt id. %")
    # draw heatmap color scale
    canvas.setFont(rFont, NfSize)
    hk_list = sorted(idpt.iterkeys(), reverse=True)
    hk_list.insert(0, 100)
    hk_i = 0
    hkY -= hk_boxX*1.5
    canvas.drawCentredString(hkX+hk_boxX+incrN*3, hkY, str(100))
    while hk_i < len(hk_list)-1:
        hk_boxY = (hk_list[hk_i]-hk_list[hk_i+1])*hk_u
        hkY -= hk_boxY
        canvas.setFillColor(HexColor(idpt[hk_list[hk_i+1]]))
        canvas.rect(hkX, hkY, hk_boxX, hk_boxY, fill=1)
        canvas.setFillColor(black)
        canvas.drawCentredString(hkX+hk_boxX+incrN*3, hkY,
                                 str(hk_list[hk_i+1]))
        hk_i += 1







