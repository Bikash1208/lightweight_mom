import ast
p = 'tune_c4.py'
s = open(p, encoding='utf-8').read()

def rep(old, new, count=1):
    global s
    assert s.count(old) == count, (s.count(old), old[:80])
    s = s.replace(old, new, 1)

# 1. loop guard
rep('v.get("gv", 0.0) < 0.5) else ()):', 'v.get("gv", 0.0) < 0.5 and v.get("bs", 0.0) < 0.5) else ()):')

# 2. BS cuts before the SLOT_SUBTRACT (only the build_c4 occurrence, which is indented 4 + 'if cuts:')
rep('''    if cuts:
        require(m.subtract(disc.name, [c.name for c in cuts], keep_originals=False), "SLOT_SUBTRACT")
''', '''    if v.get("bs", 0.0) > 0.5:
        sw_ = v["bs_strip"]; gp_ = v["bs_gap"]; L2_ = v["slot_l"] / 2; hs_ = sw_ / 2
        cuts.append(box("CUT_SA", [-w_/2, -L2_, -ct], [w_, L2_ - hs_, 3*ct], "vacuum"))
        cuts.append(box("CUT_SB", [-w_/2, hs_, -ct], [w_, L2_ - hs_, 3*ct], "vacuum"))
        cuts.append(box("CUT_GAP", [-gp_/2, -hs_, -ct], [gp_, sw_, 3*ct], "vacuum"))
        if v["bs_rim"] > 0.5:
            rc_ = v["bs_rimcut"]; ext_ = v["disc_r"] + 0.3 - L2_
            cuts.append(box("CUT_RIMP", [-rc_/2, L2_, -ct], [rc_, ext_, 3*ct], "vacuum"))
            cuts.append(box("CUT_RIMM", [-rc_/2, -L2_ - ext_, -ct], [rc_, ext_, 3*ct], "vacuum"))
    if cuts:
        require(m.subtract(disc.name, [c.name for c in cuts], keep_originals=False), "SLOT_SUBTRACT")
''')

# 3. sheets before GROUND creation inside build_c4 (the one right after the DISC_VOLUME try block)
rep('''    ground = box("GROUND", ["-period/2", "-period/2", "ground_bottom"],
                 ["period", "period", "cu_t"])
    rt = box("SUBSTRATE", ["-period/2", "-period/2", "-rt_t"],
             ["period", "period", "rt_t"], "RT5880_BAR90_TUNE")
    bias = box("BIAS_RT5880", ["-period/2", "-period/2", "bias_bottom"],
               ["period", "period", "bias_t"], "RT5880_BAR90_TUNE")
    if v.get("gv", 0.0) > 0.5:''', '''    if v.get("bs", 0.0) > 0.5:
        val(h, "Rbs", v["bs_rs"], "ohm"); val(h, "Xbs", v["bs_xs"], "ohm")
        gsh = require(m.create_rectangle("XY", [-v["bs_gap"]/2, -v["bs_strip"]/2, ct/2], [v["bs_gap"], v["bs_strip"]], name="BS_GAPSHEET", material="vacuum"), "BS_GAPSHEET")
        require(h.assign_impedance_to_sheet(gsh.name, name="BS_DIODE", resistance="Rbs", reactance="Xbs"), "BS_DIODE_BOUNDARY")
        if v["bs_rim"] > 0.5:
            val(h, "Rrc", v["bs_rrc"], "ohm"); val(h, "Xrc", v["bs_xrc"], "ohm")
            rc_ = v["bs_rimcut"]; pad_ = v["bs_cappad"]; L2_ = v["slot_l"] / 2
            ysp = L2_ + 0.35
            c1 = require(m.create_rectangle("XY", [-rc_/2, ysp, ct/2], [rc_, pad_], name="BS_CAPSHEET_P", material="vacuum"), "BS_CAP_P")
            c2 = require(m.create_rectangle("XY", [-rc_/2, -ysp - pad_, ct/2], [rc_, pad_], name="BS_CAPSHEET_M", material="vacuum"), "BS_CAP_M")
            require(h.assign_impedance_to_sheet([c1.name, c2.name], name="BS_CAP", resistance="Rrc", reactance="Xrc"), "BS_CAP_BOUNDARY")
    ground = box("GROUND", ["-period/2", "-period/2", "ground_bottom"],
                 ["period", "period", "cu_t"])
    rt = box("SUBSTRATE", ["-period/2", "-period/2", "-rt_t"],
             ["period", "period", "rt_t"], "RT5880_BAR90_TUNE")
    bias = box("BIAS_RT5880", ["-period/2", "-period/2", "bias_bottom"],
               ["period", "period", "bias_t"], "RT5880_BAR90_TUNE")
    if v.get("gv", 0.0) > 0.5:''')

# 4. cfg_bs before build_c4
rep('def build_c4(h, v, output, floor_kind="radiation"):', '''def cfg_bs(name, floor_kind):
    """Bias-network screen cell: BS_<REF|RIM|RIMX>_<ONA|ONB|OFFD>. Anchor disc r5.45, slot 8.72 (y) x 2.91 (x), 0.5 mm strip across the slot
    centre with a 0.3 mm gap carrying a 0.3(x) x 0.5(y) mm impedance sheet (Zs = Z x 0.5/0.3). RIM/RIMX: two 0.3 mm cuts across the rim at the slot
    ends (x=0, y +-4.36..5.45) each bridged by a 0.3 x 0.5 mm cap sheet (RIM: 0.3 ohm -> 0.5 ohm/sq; RIMX: +j11 ohm -> j18.33 ohm/sq)."""
    import re as _r
    m_ = _r.match(r"BS_(REF|RIM|RIMX)_(ONA|ONB|OFFD)$", name)
    require(m_ is not None, "BS_NAME")
    kind_, st_ = m_.groups()
    v = cfg_st("ST_R5p45_L8p72_W2p91", floor_kind)
    v.update(mn_n=0.0, mn_t=1.2, mn_d=0.0, mn_w=2.91, mn_edge=1.0, mn_ratio=0.0, mn_exp=2.0/8, mn_bridge=0.0, c4_s0=0.0, c4_on=0.0, c4x=0.0,
             bs=1.0, bs_strip=0.5, bs_gap=0.3, bs_rim=0.0 if kind_ == "REF" else 1.0, bs_rimcut=0.3, bs_cappad=0.5)
    rs_, xs_ = {"ONA": (3.333, -232.667), "ONB": (37.05, -559.0), "OFFD": (2.667, 48.667)}[st_]
    v.update(bs_rs=rs_, bs_xs=xs_, bs_rrc=0.5, bs_xrc=18.333 if kind_ == "RIMX" else 0.0)
    return v


def build_c4(h, v, output, floor_kind="radiation"):''')

# 5. routing
rep('or a.config.startswith("GD_")', 'or a.config.startswith("GD_") or a.config.startswith("BS_")')
rep('v=(cfg_gd(name,a.floor) if name.startswith("GD_")', 'v=(cfg_bs(name,a.floor) if name.startswith("BS_") else cfg_gd(name,a.floor) if name.startswith("GD_")')
rep('or name.startswith("GD_")) else build_mn', 'or name.startswith("GD_") or name.startswith("BS_")) else build_mn')
ast.parse(s)
open(p, 'w', encoding='utf-8').write(s)
print('patched, parses')
