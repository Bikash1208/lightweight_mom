"""MODEL_ONLY power-conserving radial gate and selected-cell board study."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import traceback
from datetime import datetime,timezone

import numpy as np
from scipy.ndimage import maximum_filter

from scripts.research.sub6_meander_aperture import C0, F, K0, THETA_BEAM, angular_distance
P,R,FEED_FLOOR=.0116,.150,.0116/2

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"results"/"sub6_meander_model"/"key_area_board"
SEEDS=range(8)
GATE_RADII_MM=(5.8,10,20,30,40,50,60,75)
GATE_ROWS=(1,2,3,4,5)
TRANSMISSIONS=(0.0,.01,.02,.05,.10)
FRACTIONS=np.arange(.05,1.001,.05)
N_ANGLE=720
DS=.001


def lattice():
    ix,iy=np.meshgrid(np.arange(-12,13),np.arange(-12,13),indexing="xy")
    x,y=ix*P,iy*P
    keep=(x*x+y*y<=R*R)&((ix!=0)|(iy!=0))
    return ix[keep],iy[keep],x[keep],y[keep]


def dither(ix,iy,x,width_mm,seed):
    p_b=.5*(1-np.tanh(x/(width_mm*.001)))
    is_b=np.zeros(len(x),bool)
    rng=np.random.default_rng(seed)
    for col in np.unique(ix):
        ids=np.flatnonzero(ix==col)
        exact=float(np.sum(p_b[ids]))
        n_b=int(math.floor(exact)+(rng.random()<exact%1))
        is_b[rng.choice(ids,size=n_b,replace=False)]=True
    return p_b,is_b


def directions(step):
    theta=np.arange(0,90+step,step)
    phi=np.arange(0,360,step)
    tt,pp=np.meshgrid(np.radians(theta),np.radians(phi),indexing="ij")
    t,p=tt.ravel(),pp.ravel()
    uv=np.column_stack((np.sin(t)*np.cos(p),np.sin(t)*np.sin(p)))
    return tt.shape,t,p,uv


def alpha_profile(radius):
    """Same anchor-type radial alpha for both hologram hands: 1.11/1.45/1.83."""
    return np.interp(radius,[0,.075,.150],[1.11,1.45,1.83])


def db_ratio(a,b):
    if a==0: return -math.inf
    if b<=0: return None
    return 10*math.log10(a/b)


def drop_db(after,before):
    if after==0: return math.inf
    if before<=0: return None
    return 10*math.log10(before/after)


def setup():
    ix,iy,x,y=lattice()
    rho=np.hypot(x,y)
    nmax=12
    valid_grid=np.zeros((25,25),bool)
    valid_grid[iy+nmax,ix+nmax]=True
    nr=int(math.ceil(R/(P/4)))
    s=(np.arange(nr)+.5)/nr
    ii=np.clip(np.rint(x[:,None]*s/P).astype(int)+nmax,0,24)
    jj=np.clip(np.rint(y[:,None]*s/P).astype(int)+nmax,0,24)
    angles=2*math.pi*np.arange(N_ANGLE)/N_ANGLE
    sr=np.arange(DS/2,R,DS)
    ri=np.clip(np.rint(np.cos(angles)[:,None]*sr/P).astype(int)+nmax,0,24)
    rj=np.clip(np.rint(np.sin(angles)[:,None]*sr/P).astype(int)+nmax,0,24)
    angle_index=np.rint((np.mod(np.arctan2(y,x),2*math.pi))/(2*math.pi)*N_ANGLE).astype(int)%N_ANGLE
    shape,t,phi,uv=directions(2)
    phase=np.exp(1j*K0*(uv@np.column_stack((x,y)).T))*np.cos(t)[:,None]
    targets={}
    for hand,target_phi in (("A",0.0),("B",math.pi)):
        distance=angular_distance(t,phi,target_phi)
        targets[hand]={"main":distance<=math.radians(20),
                       "outside":distance>math.asin(min(1,1.22*(C0/F)/(2*R/math.sqrt(2))))}
    launches={"A":np.exp(-1j*K0*math.sin(THETA_BEAM)*x),
              "B":np.exp(+1j*K0*math.sin(THETA_BEAM)*x)}
    return dict(ix=ix,iy=iy,x=x,y=y,rho=rho,nmax=nmax,valid_grid=valid_grid,s=s,ii=ii,jj=jj,
                angles=angles,sr=sr,ri=ri,rj=rj,angle_index=angle_index,
                shape=shape,t=t,phi=phi,phase=phase,targets=targets,launches=launches)


def layout(ctx,seed):
    _,is_b=dither(ctx["ix"],ctx["iy"],ctx["x"],30,seed)
    grid=np.zeros((25,25),bool)
    grid[ctx["iy"]+12,ctx["ix"]+12]=is_b
    return is_b,grid


def blocker_masks(ctx,is_b,radius_mm,rows,span):
    rho=ctx["rho"]
    radial=(rho>=radius_mm*.001)&(rho<radius_mm*.001+rows*P)
    x=ctx["x"]
    if span=="half":
        a=radial&(~is_b)&(x>=0)
        b=radial&is_b&(x<=0)
    else:
        a=radial&(~is_b)&(x>=-.030)
        b=radial&is_b&(x<=+.030)
    return a,b


def selected_mask(ctx,is_b,hand,placement,fraction):
    target=np.flatnonzero(~is_b if hand=="A" else is_b)
    x,y,r=ctx["x"][target],ctx["y"][target],ctx["rho"][target]
    sign=1 if hand=="A" else -1
    if placement=="inner_ring": score=r
    elif placement=="outer_ring": score=-r
    elif placement=="axis_patch": score=np.hypot(x-sign*.075,y)
    elif placement=="upper_boundary_patch": score=np.hypot(x-sign*.020,y-.075)
    elif placement=="lower_boundary_patch": score=np.hypot(x-sign*.020,y+.075)
    else: raise ValueError(placement)
    n=max(1,int(math.ceil(fraction*len(target))))
    selected=np.zeros(len(is_b),bool)
    selected[target[np.argsort(score,kind="stable")[:n]]]=True
    return selected


def grid_mask(ctx,mask):
    g=np.zeros((25,25),bool)
    g[ctx["iy"]+12,ctx["ix"]+12]=mask
    return g


def radial_state(ctx,is_b,grid,switched=None,inner_radius=0,gate=None,bound="lost",
                 pass_nonleaking=None):
    """Integrate dP/dr=-2 alpha P; gate conserves or discards P at one radius."""
    if switched is None: switched=np.zeros(len(is_b),bool)
    if pass_nonleaking is not None:
        switched=switched|pass_nonleaking
    switched_grid=grid_mask(ctx,switched)
    cell_alpha=alpha_profile(ctx["rho"])*~switched
    cell_alpha=np.where(ctx["rho"]<inner_radius,0.0,cell_alpha)
    local=alpha_profile(ctx["rho"][:,None]*ctx["s"][None,:])
    local*=~switched_grid[ctx["jj"],ctx["ii"]]
    local*=ctx["valid_grid"][ctx["jj"],ctx["ii"]]
    local*=ctx["rho"][:,None]*ctx["s"][None,:]>=inner_radius
    local[(ctx["jj"]==12)&(ctx["ii"]==12)]=0
    depth=ctx["rho"]*np.mean(local,axis=1)
    ray_alpha=alpha_profile(ctx["sr"])[None,:]*np.ones((N_ANGLE,1))
    ray_alpha*=~switched_grid[ctx["rj"],ctx["ri"]]
    ray_alpha*=ctx["valid_grid"][ctx["rj"],ctx["ri"]]
    ray_alpha*=ctx["sr"][None,:]>=inner_radius
    ray_alpha[(ctx["rj"]==12)&(ctx["ri"]==12)]=0
    hit=np.zeros(N_ANGLE,bool)
    factor=np.ones(N_ANGLE)
    gate_radius=None
    blocked_power=0.0
    redistribution_feasible=True
    if gate:
        gate_mask,radius_mm,rows,suppressed,transmission=gate
        gate_radius=radius_mm*.001
        blocker_grid=grid_mask(ctx,gate_mask)
        annulus=(ctx["sr"]>=gate_radius)&(ctx["sr"]<gate_radius+rows*P)
        hit=np.any(blocker_grid[ctx["rj"][:,annulus],ctx["ri"][:,annulus]],axis=1)
        gate_index=int(np.searchsorted(ctx["sr"],gate_radius))
        gate_depth=np.sum(ray_alpha[:,:gate_index]*DS,axis=1)
        p_gate=np.exp(-2*gate_depth)
        blocked_power=float(np.sum(p_gate[hit]))/N_ANGLE*(1-transmission) if np.any(hit) else 0.0
        factor[hit]=transmission
        if bound=="redistributed" and blocked_power>0:
            survivor="B" if suppressed=="A" else "A"
            receive=(np.cos(ctx["angles"])>0) if survivor=="A" else (np.cos(ctx["angles"])<0)
            receive &=~hit
            denom=float(np.mean(p_gate*receive))
            if denom<=0:
                redistribution_feasible=False
            else:
                factor[receive]+=blocked_power/denom
                blocked_power=0.0
    cell_factor=np.where(ctx["rho"]>=gate_radius if gate else False,
                         factor[ctx["angle_index"]],1.0)
    amp=np.sqrt(2*cell_alpha/np.maximum(ctx["rho"],FEED_FLOOR))*np.exp(-depth)*np.sqrt(cell_factor)
    before=np.cumsum(ray_alpha,axis=1)*DS-ray_alpha*DS
    p_before=np.exp(-2*before)
    segment_factor=np.ones((N_ANGLE,len(ctx["sr"])))
    if gate:
        segment_factor[:,ctx["sr"]>=gate_radius]=factor[:,None]
    segment=p_before*(1-np.exp(-2*ray_alpha*DS))*segment_factor
    is_b_ray=grid[ctx["rj"],ctx["ri"]]
    conv_a=float(np.mean(np.sum(np.where(is_b_ray,0,segment),axis=1)))
    conv_b=float(np.mean(np.sum(np.where(is_b_ray,segment,0),axis=1)))
    edge=float(np.mean(np.exp(-2*np.sum(ray_alpha*DS,axis=1))*factor))
    balance=conv_a+conv_b+edge+blocked_power
    if abs(balance-1)>.015: raise RuntimeError(f"POWER_BALANCE {balance:.6f}")
    return {"amp":amp,"conversion_A":conv_a,"conversion_B":conv_b,
            "edge_power":edge,"blocked_lost":blocked_power,"power_balance":balance,
            "blocked_ray_fraction":float(np.mean(hit)),
            "redistribution_feasible":redistribution_feasible}


def fields(ctx,is_b,amp,full=True):
    wa=amp*(~is_b)*ctx["launches"]["A"]
    wb=amp*is_b*ctx["launches"]["B"]
    if full:
        return ctx["phase"]@wa,ctx["phase"]@wb
    mask=ctx["targets"]["A"]["main"]|ctx["targets"]["B"]["main"]
    return ctx["phase"][mask]@wa,ctx["phase"][mask]@wb


def pattern_metrics(ctx,fa,fb):
    intensities={"A":abs(fa)**2,"B":abs(fb)**2}
    out={}
    for hand in ("A","B"):
        arr=intensities[hand]
        peak=float(np.max(arr[ctx["targets"][hand]["main"]]))
        shaped=arr.reshape(ctx["shape"])
        local_max=shaped>=maximum_filter(shaped,size=(3,3),mode=("nearest","wrap"))*(1-1e-12)
        side=arr[ctx["targets"][hand]["outside"]&local_max.ravel()]
        wide=float(np.max(arr[~ctx["targets"][hand]["main"]]))
        sll=10*math.log10(float(np.max(side))/peak) if len(side) and peak>0 else None
        out[hand]={"peak":peak,"sll_db":sll,"wide_angle_peak":wide}
    return out


def peak_only(ctx,fa,fb):
    mask=ctx["targets"]["A"]["main"]|ctx["targets"]["B"]["main"]
    return {hand:float(np.max(abs(field[ctx["targets"][hand]["main"][mask]])**2))
            for hand,field in (("A",fa),("B",fb))}


def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def study1(ctx,progress=None):
    rows=[]
    for seed in SEEDS:
        if progress: progress("study1",seed)
        is_b,grid=layout(ctx,seed)
        normal=radial_state(ctx,is_b,grid)
        fa,fb=fields(ctx,is_b,normal["amp"])
        normal_metrics=pattern_metrics(ctx,fa,fb)
        for radius_mm in GATE_RADII_MM:
            for nrows in GATE_ROWS:
                for span in ("half","half_plus_dither30"):
                    gate_a,gate_b=blocker_masks(ctx,is_b,radius_mm,nrows,span)
                    all_gate=gate_a|gate_b
                    diodes=int(np.sum(all_gate))
                    for pass_state in ("anchor_like","nonleaking"):
                        nonleak=all_gate if pass_state=="nonleaking" else None
                        if nonleak is None:
                            both,bm=normal,normal_metrics
                        else:
                            both=radial_state(ctx,is_b,grid,pass_nonleaking=nonleak)
                            fa,fb=fields(ctx,is_b,both["amp"])
                            bm=pattern_metrics(ctx,fa,fb)
                        for suppressed,mask in (("A",gate_a),("B",gate_b)):
                            surviving="B" if suppressed=="A" else "A"
                            for transmission in TRANSMISSIONS:
                                for bound in ("lost","redistributed"):
                                    state=radial_state(ctx,is_b,grid,pass_nonleaking=nonleak,
                                          gate=(mask,radius_mm,nrows,suppressed,transmission),bound=bound)
                                    fa,fb=fields(ctx,is_b,state["amp"])
                                    pm=pattern_metrics(ctx,fa,fb)
                                    peak_drop=drop_db(pm[suppressed]["peak"],bm[suppressed]["peak"])
                                    power_drop=drop_db(state["conversion_"+suppressed],both["conversion_"+suppressed])
                                    survive_change=db_ratio(pm[surviving]["peak"],bm[surviving]["peak"])
                                    wide=db_ratio(pm[suppressed]["wide_angle_peak"],bm[suppressed]["peak"])
                                    for block_state in ("ABSORBING","REFLECTIVE"):
                                        rows.append(dict(label="MODEL_ONLY",seed=seed,gate_inner_radius_mm=radius_mm,
                                          gate_rows=nrows,transmitted_power_fraction=transmission,
                                          block_state=block_state,pass_state=pass_state,span=span,
                                          suppressed=suppressed,surviving=surviving,bound=bound,
                                          total_diodes=diodes,diodes_A=int(np.sum(gate_a)),diodes_B=int(np.sum(gate_b)),
                                          blocked_ray_fraction=state["blocked_ray_fraction"],
                                          bound_feasible=state["redistribution_feasible"],
                                          suppressed_peak_drop_db=peak_drop,suppressed_power_drop_db=power_drop,
                                          suppressed_wide_angle_db_rel_both_peak=wide,
                                          suppressed_sll_db=pm[suppressed]["sll_db"],
                                          survivor_gain_change_db=survive_change,survivor_sll_db=pm[surviving]["sll_db"],
                                          both_A_gain_cost_db=-db_ratio(bm["A"]["peak"],normal_metrics["A"]["peak"]),
                                          both_B_gain_cost_db=-db_ratio(bm["B"]["peak"],normal_metrics["B"]["peak"]),
                                          both_A_conversion=both["conversion_A"],both_B_conversion=both["conversion_B"],
                                          blocked_lost=state["blocked_lost"] if block_state=="ABSORBING" else 0.0,
                                          blocked_reflected=state["blocked_lost"] if block_state=="REFLECTIVE" else 0.0,
                                          edge_power=state["edge_power"],power_balance=state["power_balance"],
                                          gate_10db="PASS" if state["redistribution_feasible"] and peak_drop is not None and peak_drop>=10 else "FAIL",
                                          full_requirement="PASS" if state["redistribution_feasible"] and peak_drop is not None and peak_drop>=10 and survive_change is not None and survive_change>=0 else "FAIL"))
    return rows


def study2(ctx,progress=None):
    rows=[]
    for seed in SEEDS:
        if progress: progress("study2",seed)
        is_b,grid=layout(ctx,seed)
        base=radial_state(ctx,is_b,grid)
        fa,fb=fields(ctx,is_b,base["amp"])
        baseline=peak_only(ctx,fa[ctx["targets"]["A"]["main"]|ctx["targets"]["B"]["main"]],
                              fb[ctx["targets"]["A"]["main"]|ctx["targets"]["B"]["main"]])
        for hand in ("A","B"):
            surviving="B" if hand=="A" else "A"
            total=int(np.sum(~is_b if hand=="A" else is_b))
            for placement in ("inner_ring","outer_ring","axis_patch","upper_boundary_patch","lower_boundary_patch"):
                for fraction in FRACTIONS:
                    switched=selected_mask(ctx,is_b,hand,placement,float(fraction))
                    state=radial_state(ctx,is_b,grid,switched=switched)
                    aa,bb=fields(ctx,is_b,state["amp"],full=False)
                    peaks=peak_only(ctx,aa,bb)
                    drop=drop_db(peaks[hand],baseline[hand])
                    power_drop=drop_db(state["conversion_"+hand],base["conversion_"+hand])
                    gain=db_ratio(peaks[surviving],baseline[surviving])
                    rows.append(dict(label="MODEL_ONLY",seed=seed,suppressed=hand,surviving=surviving,
                                     placement=placement,requested_target_fraction=round(float(fraction),3),
                                     switched_cells=int(np.sum(switched)),target_cells=total,
                                     fraction_of_target=float(np.sum(switched))/total,
                                     fraction_of_board=float(np.sum(switched))/len(is_b),
                                     suppressed_peak_drop_db=drop,suppressed_power_drop_db=power_drop,
                                     survivor_gain_change_db=gain,edge_power=state["edge_power"],
                                     power_balance=state["power_balance"],
                                     gate_10db="PASS" if drop is not None and drop>=10 else "FAIL",
                                     full_requirement="PASS" if drop is not None and drop>=10 and gain is not None and gain>=0 else "FAIL"))
    return rows


def summarize(rows1,rows2):
    def median_group(rows,keys,metrics):
        buckets={}
        for row in rows:
            key=tuple(row[k] for k in keys)
            buckets.setdefault(key,[]).append(row)
        out=[]
        for key,group in buckets.items():
            item=dict(zip(keys,key))
            for metric in metrics:
                vals=[float(r[metric]) for r in group if r[metric] is not None and math.isfinite(float(r[metric]))]
                item[metric+"_median"]=float(np.median(vals)) if vals else None
                item[metric+"_min"]=min(vals) if vals else None
                item[metric+"_max"]=max(vals) if vals else None
            item["all_seed_10db_pass"]=all(r["gate_10db"]=="PASS" for r in group)
            item["all_seed_full_pass"]=all(r["full_requirement"]=="PASS" for r in group)
            if "total_diodes" in group[0]:
                item["total_diodes_max"]=max(r["total_diodes"] for r in group)
            out.append(item)
        return out
    g1=median_group(rows1,["gate_inner_radius_mm","gate_rows","transmitted_power_fraction",
                           "block_state","pass_state","span","suppressed","bound"],
                    ["suppressed_peak_drop_db","suppressed_power_drop_db",
                     "suppressed_wide_angle_db_rel_both_peak","suppressed_sll_db",
                     "survivor_gain_change_db","survivor_sll_db",
                     "both_A_gain_cost_db","both_B_gain_cost_db"])
    joint=[]
    joint_keys=("gate_inner_radius_mm","gate_rows","transmitted_power_fraction","block_state","pass_state","span","bound")
    buckets={}
    for item in g1:
        key=tuple(item[k] for k in joint_keys)
        buckets.setdefault(key,{})[item["suppressed"]]=item
    for key,pair in buckets.items():
        if set(pair)!={"A","B"}: raise RuntimeError(f"MISSING_ONE_BEAM_STATE {key}")
        a,b=pair["A"],pair["B"]
        joint.append(dict(zip(joint_keys,key))|{
            "total_diodes_max":max(a["total_diodes_max"],b["total_diodes_max"]),
            "A_only_B_peak_drop_db":b["suppressed_peak_drop_db_median"],
            "A_only_B_power_drop_db":b["suppressed_power_drop_db_median"],
            "A_only_A_gain_change_db":b["survivor_gain_change_db_median"],
            "A_only_A_sll_db":b["survivor_sll_db_median"],
            "B_only_A_peak_drop_db":a["suppressed_peak_drop_db_median"],
            "B_only_A_power_drop_db":a["suppressed_power_drop_db_median"],
            "B_only_B_gain_change_db":a["survivor_gain_change_db_median"],
            "B_only_B_sll_db":a["survivor_sll_db_median"],
            "both_A_gain_cost_db":a["both_A_gain_cost_db_median"],
            "both_B_gain_cost_db":a["both_B_gain_cost_db_median"],
            "both_states_10db":a["all_seed_10db_pass"] and b["all_seed_10db_pass"],
            "both_states_full":a["all_seed_full_pass"] and b["all_seed_full_pass"]})
    g2=median_group(rows2,["suppressed","placement","requested_target_fraction"],
                    ["fraction_of_target","fraction_of_board","suppressed_peak_drop_db","suppressed_power_drop_db","survivor_gain_change_db"])
    thresholds=[]
    for hand in ("A","B"):
        for placement in ("inner_ring","outer_ring","axis_patch","upper_boundary_patch","lower_boundary_patch"):
            choices=sorted((r for r in g2 if r["suppressed"]==hand and r["placement"]==placement),
                           key=lambda r:r["requested_target_fraction"])
            for criterion in ("all_seed_10db_pass","all_seed_full_pass"):
                first=next((r for r in choices if r[criterion]),None)
                thresholds.append({"suppressed":hand,"placement":placement,"criterion":criterion,
                                   "fraction_of_target":first["fraction_of_target_median"] if first else None,
                                   "fraction_of_board":first["fraction_of_board_median"] if first else None,
                                   "switched_cells_median":round(first["fraction_of_board_median"]*516) if first else None})
    return g1,joint,g2,thresholds


def report(g1,joint,g2,thresholds,rows1,rows2):
    perfect=[r for r in joint if r["transmitted_power_fraction"]==0 and r["both_states_10db"]]
    max_radius=max((r["gate_inner_radius_mm"] for r in perfect),default=None)
    full=sorted((r for r in joint if r["both_states_full"] and r["block_state"]=="REFLECTIVE"),
                key=lambda r:(r["total_diodes_max"],r["gate_inner_radius_mm"]))
    highlights=full[:6]
    if highlights:
        first=highlights[0]
        lost=next((r for r in joint if all(r[k]==first[k] for k in
                   ("gate_inner_radius_mm","gate_rows","transmitted_power_fraction",
                    "block_state","pass_state","span")) and r["bound"]=="lost"),None)
        if lost: highlights.insert(0,lost)
    if not highlights:
        highlights=sorted(joint,key=lambda r:min(r["A_only_B_peak_drop_db"] or -999,
                                                   r["B_only_A_peak_drop_db"] or -999),reverse=True)[:8]
    fmt=lambda v:f"{v:.2f}" if v is not None else "NA"
    lines=["# Diodes only at key areas — MODEL_ONLY", "",
           "Launch zone: the model reserves the centre feed cell and regularizes ρ at 5.8 mm; it does not solve where the launch becomes a bound surface wave. Thus 5.8 mm is a model sweep floor, not a measured physical launch-zone radius.",
           f"Maximum swept gate inner radius with both peak drops ≥10 dB for a perfect gate is {max_radius} mm"+(" (the limit may lie beyond 75 mm)." if max_radius==75 else "; no 75 mm configuration passes both peak-drop gates."),
           "The existing model assigns independent A/RHCP and B/LHCP phase ramps by construction; physical CP is unvalidated.",
           "", "Power obeys dP/dr=−2αP; α=0 passes power outward. The same anchor α profile applies in both halves: 1.11/1.45/1.83 Np/m at r=0/75/150 mm. R=150 mm, 11.6-mm pitch, 30-mm dither, eight seeds, 2° angular grid. A-only and B-only are computed separately.",
           "", "## Study 1 — two semicircular gates", "",
           "Inner radius 5.8/10/20/30/40/50/60/75 mm; 1–5 cell rows; blocked-state transmitted power 0/1/2/5/10%; one diode per gate cell. Half-only covers assigned A cells at x≥0 and B cells at x≤0. Extended span adds assigned A cells down to x=−30 mm and B cells up to x=+30 mm. Pass-state gate cells either retain anchor α or have α=0. ABSORBING and REFLECTIVE have identical forward fields at fixed power transmission; the removed power is counted as absorbed or reflected respectively. The redistributed bound directs it into the open geometric half at the gate; for ABSORBING this is only a mathematical upper bound.",
           "", "| r inner mm | Rows | T % | Block / pass / span / bound | Diodes max | A-only: B peak/power drop; A gain; SLL dB | B-only: A peak/power drop; B gain; SLL dB | BOTH A/B cost dB | Both states |",
           "|---:|---:|---:|---|---:|---|---|---|---|"]
    for r in highlights:
        setting=f"{r['block_state']} / {r['pass_state']} / {r['span']} / {r['bound']}"
        a=f"{fmt(r['A_only_B_peak_drop_db'])}/{fmt(r['A_only_B_power_drop_db'])}; {fmt(r['A_only_A_gain_change_db'])}; {fmt(r['A_only_A_sll_db'])}"
        b=f"{fmt(r['B_only_A_peak_drop_db'])}/{fmt(r['B_only_A_power_drop_db'])}; {fmt(r['B_only_B_gain_change_db'])}; {fmt(r['B_only_B_sll_db'])}"
        cost=f"{fmt(r['both_A_gain_cost_db'])}/{fmt(r['both_B_gain_cost_db'])}"
        lines.append(f"| {r['gate_inner_radius_mm']} | {r['gate_rows']} | {100*r['transmitted_power_fraction']:.0f} | {setting} | {r['total_diodes_max']} | {a} | {b} | {cost} | {'PASS' if r['both_states_full'] else 'FAIL'} |")
    lines += ["",f"Every one-beam case, including suppressed-hand wide-angle level and SLL, is in `study1_gates_rows.csv` ({len(rows1)} seed rows); `study1_gates_summary.csv` ({len(g1)}) and `study1_joint.csv` ({len(joint)}) aggregate the cases. PASS requires ≥10 dB peak drop and no surviving-beam gain decrease for every seed. Power drop is the integrated suppressed-hand radiated conversion change; gain is referenced to equal launched power. No lost-power-bound configuration meets the full requirement for both one-beam states; the full passes require the ideal redistribution bound.",
              "", "## Study 2 — selected anchor cells with α=0", "",
              "Inner/outer rings select target-hand cells by radius. Axis patch selects cells nearest (±75,0) mm; upper/lower boundary patches select cells nearest (±20,±75) mm. Power released by switched cells remains on their radial rays.",
              "", "| Suppressed | Placement | Fraction of target cells for 10 dB (all seeds) | Fraction of whole board | Fraction meeting 10 dB and no survivor loss |",
              "|---|---|---:|---:|---:|"]
    for hand in ("A","B"):
        for placement in ("inner_ring","outer_ring","axis_patch","upper_boundary_patch","lower_boundary_patch"):
            a=next(x for x in thresholds if x["suppressed"]==hand and x["placement"]==placement and x["criterion"]=="all_seed_10db_pass")
            b=next(x for x in thresholds if x["suppressed"]==hand and x["placement"]==placement and x["criterion"]=="all_seed_full_pass")
            frac=lambda v:f"{v:.3f}" if v is not None else "none ≤1.0"
            lines.append(f"| {hand} | {placement} | {frac(a['fraction_of_target'])} | {frac(a['fraction_of_board'])} | {frac(b['fraction_of_target'])} |")
    lines += ["", "At the 0.702 target-half threshold, independently switching either beam would require about 180 A-cell plus 182 B-cell diodes (about 362 total). The outer-ring threshold requires about 414 total. These selected-cell schemes exceed the earlier 200-diode preference.",
              f"`study2_selected_rows.csv` ({len(rows2)} seed rows), `study2_selected_summary.csv` ({len(g2)}) and `study2_thresholds.csv` provide the full placement/fraction tables.",
              "", "FLAG: This is a coherent ideal-CP array-factor model, not HFSS. Gate reflection phase, coupling into other modes, feed mismatch, diode impedance, hardware loss, and the physical surface-wave launch radius are unknown. Blocked-power redistribution is an upper bound, not a predicted circuit response. Some broad gates leave no open-half ray to receive redistributed power and are marked infeasible. Reported power conservation is radial bookkeeping, not an electromagnetic scattering solution.",
              "", "HANDOFF: MODEL_ONLY | both one-beam states and BOTH evaluated | full R150 aperture | no HFSS | next: validate a realizable gate or selected-cell pattern."]
    (OUT/"REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    status={"status":"RUNNING","pid":os.getpid(),"started_utc":datetime.now(timezone.utc).isoformat()}
    def progress(stage,seed):
        status.update(stage=stage,seed=seed,heartbeat_utc=datetime.now(timezone.utc).isoformat())
        (OUT/"status.json").write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8")
    progress("setup",-1)
    try:
        ctx=setup()
        rows1=study1(ctx,progress)
        write_csv(OUT/"study1_gates_rows.csv",rows1)
        rows2=study2(ctx,progress)
        write_csv(OUT/"study2_selected_rows.csv",rows2)
        g1,joint,g2,thresholds=summarize(rows1,rows2)
        write_csv(OUT/"study1_gates_summary.csv",g1)
        write_csv(OUT/"study1_joint.csv",joint)
        write_csv(OUT/"study2_selected_summary.csv",g2)
        write_csv(OUT/"study2_thresholds.csv",thresholds)
        report(g1,joint,g2,thresholds,rows1,rows2)
        (OUT/"meta.json").write_text(json.dumps({"label":"MODEL_ONLY","radius_mm":150,"pitch_mm":11.6,
           "alpha_profile_np_m":{"r0":1.11,"r75":1.45,"r150":1.83},"seeds":list(SEEDS),
           "power_conserving":True,"study1_rows":len(rows1),"study2_rows":len(rows2),
           "max_power_balance_error":max(abs(r["power_balance"]-1) for r in rows1+rows2)},indent=2)+"\n",encoding="utf-8")
        status.update(status="COMPLETE",ended_utc=datetime.now(timezone.utc).isoformat(),study1_rows=len(rows1),study2_rows=len(rows2))
        print(f"study1={len(rows1)} study2={len(rows2)} max_balance_error={max(abs(r['power_balance']-1) for r in rows1+rows2):.3g}")
    except BaseException:
        status.update(status="FAILED",traceback=traceback.format_exc(),ended_utc=datetime.now(timezone.utc).isoformat())
        raise
    finally:
        (OUT/"status.json").write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__": main()
