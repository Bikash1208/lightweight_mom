import ast
p = 'tune_c4.py'
s = open(p, encoding='utf-8').read()

def rep(old, new):
    global s
    assert s.count(old) == 1, (s.count(old), old[:90])
    s = s.replace(old, new, 1)

# cfg_bs: add VIA kind with choke type
rep('m_ = _r.match(r"BS_(REF|RIM|RIMX)_(ONA|ONB|OFFD)$", name)\n    require(m_ is not None, "BS_NAME")\n    kind_, st_ = m_.groups()',
    'm_ = _r.match(r"BS_(REF|RIM|RIMX|VIA)_(ONA|ONB|OFFD)(?:_(J1000|R1000))?$", name)\n    require(m_ is not None, "BS_NAME")\n    kind_, st_, ch_ = m_.groups()')
rep('bs_rim=0.0 if kind_ == "REF" else 1.0, bs_rimcut=0.3, bs_cappad=0.5)',
    'bs_rim=0.0 if kind_ == "REF" else 1.0, bs_rimcut=0.3, bs_cappad=0.5, bs_via=1.0 if kind_ == "VIA" else 0.0,\n             bs_rv=1000.0 if ch_ == "R1000" else 1.0, bs_xv=0.0 if ch_ == "R1000" else 1000.0)')

via = '''
    if v.get("bs", 0.0) > 0.5 and v.get("bs_via", 0.0) > 0.5:
        rt_t_ = v["rt_t"]; z0 = -rt_t_ - v["cu_t"]; zt, zb = -1.6, -1.9
        val(h, "Rvb", v["bs_rv"], "ohm"); val(h, "Xvb", v["bs_xv"], "ohm")
        shs = []
        for tag, gx in (("P", 3.3), ("M", -3.3)):
            hole = box("BV_HOLE_" + tag, [gx - .15, -.15, z0], [.3, .3, rt_t_ + v["cu_t"]], "vacuum")
            require(m.subtract(rt.name, hole.name, keep_originals=False), "BV_HOLE_SUB_" + tag)
            up = box("BV_UP_" + tag, [gx - .15, -.15, zt], [.3, .3, -zt + v["cu_t"]/2])
            lo = box("BV_LO_" + tag, [gx - .15, -.15, z0], [.3, .3, zb - z0])
            lumps = [n_ for n_ in m.get_bodynames_from_position([gx, 0, v["cu_t"]/2]) if n_.startswith("DISC")]
            print("BV_LUMP", tag, lumps, flush=True)
            require(len(lumps) == 1, "BV_LUMP_" + tag)
            require(m.unite([lumps[0], up.name]), "BV_UNITE_UP_" + tag)
            require(m.unite([ground.name, lo.name]), "BV_UNITE_LO_" + tag)
            sh_ = require(m.create_rectangle("XZ", [gx - .15, 0, zb], [.3, zt - zb], name="BV_SHEET_" + tag, material="vacuum"), "BV_SHEET_" + tag)
            shs.append(sh_.name)
        require(h.assign_impedance_to_sheet(shs, name="BV_CHOKE", resistance="Rvb", reactance="Xvb"), "BV_BOUNDARY")
        _lumps = [o for o in m.object_names if o.startswith("DISC")]
        print("DISC_LUMPS_AFTER_VIAS", _lumps, flush=True)
'''
anchor = 'require(h.assign_impedance_to_sheet(sh.name, name="GV_DIODE", resistance="Rgv", reactance="Xgv"), "GV_BOUNDARY")\n'
rep(anchor, anchor + via)
ast.parse(s)
open(p, 'w', encoding='utf-8').write(s)
print('patched, parses')
