import tkinter as tk
from tkinter import font as tkfont
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


TH = {
    "bg":       "#05070f",
    "bg2":      "#0a0d1e",
    "card":     "#0e1128",
    "card2":    "#141836",
    "input":    "#1a1e40",
    "border":   "#252b5a",
    "hi":       "#7c3aed",
    "text":     "#f0f3ff",
    "text2":    "#9aa3cc",
    "text3":    "#525d94",
    "purple":   "#a855f7",
    "blue":     "#3b82f6",
    "cyan":     "#06b6d4",
    "teal":     "#14b8a6",
    "green":    "#22c55e",
    "lime":     "#84cc16",
    "yellow":   "#eab308",
    "amber":    "#f59e0b",
    "orange":   "#f97316",
    "rose":     "#f43f5e",
    "pink":     "#ec4899",
    "j1":       "#f472b6",
    "j2":       "#a78bfa",
    "j3":       "#34d399",
    "j4":       "#fbbf24",
    "link1_a":  "#6366f1",
    "link1_b":  "#818cf8",
    "link2_a":  "#0891b2",
    "link2_b":  "#22d3ee",
    "base_c":   "#4f46e5",
    "tcp_c":    "#f43f5e",
    "grip_c":   "#f59e0b",
    "grid_c":   "#12163a",
    "floor":    "#0a0e22",
    "trough":   "#1e2354",
   
    "ped_top":  "#2a2f6e",
    "ped_side": "#1e2258",
    "ped_edge": "#4a51a8",
    "ped_ring": "#6366f1",
    "mount_c":  "#3730a3",
    "mount_hi": "#818cf8",
}


BASE_H   = 80   
BASE_R   = 55   
NECK_R   = 28   
PLATE_H  = 12   


# DH Kinematics

def dh_matrix(theta_deg, d, a, alpha_deg):
    t, al = np.radians(theta_deg), np.radians(alpha_deg)
    ct, st, ca, sa = np.cos(t), np.sin(t), np.cos(al), np.sin(al)
    return np.array([
        [ct, -st*ca,  st*sa, a*ct],
        [st,  ct*ca, -ct*sa, a*st],
        [0,   sa,     ca,    d   ],
        [0,   0,      0,     1   ]])

def forward_kinematics(t1, t2, d3, t4, l1, l2):
    T01 = dh_matrix(t1,0,l1,0)
    T12 = dh_matrix(t2,0,l2,180)
    T23 = dh_matrix(0,d3,0,0)
    T34 = dh_matrix(t4,0,0,0)
    T02=T01@T12; T03=T02@T23; T04=T03@T34
    phi = t1+t2+t4; pr = np.radians(phi)
    p = [np.zeros(3), T01[:3,3], T02[:3,3], T03[:3,3], T04[:3,3]]
    tcp = p[4]

    
    dx, dy = np.cos(pr), np.sin(pr)
    
    px, py = -np.sin(pr), np.cos(pr)

    GL = 10  
    GW = 20   
    GH = 12   

    
    lb_base = [tcp[0]+px*GW, tcp[1]+py*GW, tcp[2]]
    rb_base = [tcp[0]-px*GW, tcp[1]-py*GW, tcp[2]]

    
    lb_tip = [tcp[0]+px*GW, tcp[1]+py*GW, tcp[2]-GL]
    rb_tip = [tcp[0]-px*GW, tcp[1]-py*GW, tcp[2]-GL]

    
    grip = {
        "lb":    lb_base,
        "lt":    lb_tip,
        "rb":    rb_base,
        "rt":    rb_tip,
        
        "cross_l": lb_base,
        "cross_r": rb_base,
        
        "ar": [tcp[0]+dx*30, tcp[1]+dy*30, tcp[2]],
        
        "body_top": [tcp[0], tcp[1], tcp[2]+GH],
        "body_bot": [tcp[0], tcp[1], tcp[2]],
    }
    return {"T01":T01,"T02":T02,"T03":T03,"T04":T04,"T":T04,
            "pos":p, "x":T04[0,3],"y":T04[1,3],"z":T04[2,3],
            "phi":phi,"pr":pr,"grip":grip}

def inverse_kinematics(x, y, z, phi, l1, l2):
    r2=x**2+y**2; r=np.sqrt(r2); rm=l1+l2; rn=abs(l1-l2)
    if r>rm+1e-6: return {"error":f"Ngoài tầm với! r={r:.1f}>{rm:.1f}"}
    if r<rn-1e-6 and rn>1e-6: return {"error":f"Vùng chết! r={r:.1f}<{rn:.1f}"}
    d3=-z
    if d3<-1e-6: return {"error":f"d₃={d3:.1f}<0!"}
    c2=np.clip((r2-l1**2-l2**2)/(2*l1*l2),-1,1)
    sols={}
    for n,s in [("elbow_up",1),("elbow_down",-1)]:
        th2=s*np.arccos(c2)
        th1=np.arctan2(y,x)-np.arctan2(l2*np.sin(th2),l1+l2*np.cos(th2))
        th4=np.radians(phi)-th1-th2
        t1d=((np.degrees(th1)+180)%360)-180
        t2d=((np.degrees(th2)+180)%360)-180
        t4d=((np.degrees(th4)+180)%360)-180
        fk=forward_kinematics(t1d,t2d,d3,t4d,l1,l2)
        sols[n]={"t1":round(t1d,2),"t2":round(t2d,2),"d3":round(d3,2),"t4":round(t4d,2),"T":fk["T"]}
    return sols


# TRAJECTORY PLANNING

def trapezoidal_profile(T_total, N, accel_ratio=0.25):
    """Trapezoidal velocity profile → s(t) in [0,1]."""
    t = np.linspace(0, T_total, N)
    ta = T_total * accel_ratio
    td = T_total * accel_ratio
    tc = T_total - ta - td
    v_max = 1.0 / (tc + 0.5*ta + 0.5*td)
    s = np.zeros(N)
    for i, ti in enumerate(t):
        if ti <= ta:
            s[i] = 0.5 * v_max / ta * ti**2
        elif ti <= ta + tc:
            s[i] = 0.5 * v_max * ta + v_max * (ti - ta)
        else:
            dt = ti - (ta + tc)
            s[i] = 0.5*v_max*ta + v_max*tc + v_max*dt - 0.5*v_max/td*dt**2
    s = s / s[-1]
    return t, s

def gen_line_trajectory(p_start, p_end, N=100, T=3.0):
    """Straight line in task space."""
    _, s = trapezoidal_profile(T, N)
    ps, pe = np.array(p_start), np.array(p_end)
    return np.array([ps + si*(pe - ps) for si in s])

def gen_circle_trajectory(center, radius, z=0, N=200, T=5.0):
    """Full circle in XY plane."""
    _, s = trapezoidal_profile(T, N)
    angles = 2 * np.pi * s
    pts = np.zeros((N, 3))
    pts[:, 0] = center[0] + radius * np.cos(angles)
    pts[:, 1] = center[1] + radius * np.sin(angles)
    pts[:, 2] = z
    return pts



def draw_cylinder_3d(ax, p1, p2, radius=12, color="#818cf8", alpha=0.85, n=12):
    """Draw a 3D cylinder between two points."""
    p1,p2 = np.array(p1),np.array(p2)
    v = p2-p1; length = np.linalg.norm(v)
    if length < 0.1: return
    v = v/length
    if abs(v[2]) < 0.99: not_v = np.array([0,0,1])
    else: not_v = np.array([1,0,0])
    n1 = np.cross(v, not_v); n1/=np.linalg.norm(n1)
    n2 = np.cross(v, n1)
    theta = np.linspace(0, 2*np.pi, n+1)
    t_line = np.array([0, 1])
    theta_grid, t_grid = np.meshgrid(theta, t_line)
    X = p1[0]+v[0]*t_grid*length + radius*(n1[0]*np.cos(theta_grid)+n2[0]*np.sin(theta_grid))
    Y = p1[1]+v[1]*t_grid*length + radius*(n1[1]*np.cos(theta_grid)+n2[1]*np.sin(theta_grid))
    Z = p1[2]+v[2]*t_grid*length + radius*(n1[2]*np.cos(theta_grid)+n2[2]*np.sin(theta_grid))
    ax.plot_surface(X, Y, Z, color=color, alpha=alpha, shade=True,
                    lightsource=matplotlib.colors.LightSource(azdeg=315,altdeg=45))

def draw_sphere_3d(ax, center, radius=15, color="#f472b6", alpha=0.9, n=10):
    """Draw a 3D sphere at a point."""
    u,v = np.meshgrid(np.linspace(0,2*np.pi,n), np.linspace(0,np.pi,n))
    X = center[0]+radius*np.cos(u)*np.sin(v)
    Y = center[1]+radius*np.sin(u)*np.sin(v)
    Z = center[2]+radius*np.cos(v)
    ax.plot_surface(X,Y,Z, color=color, alpha=alpha, shade=True)

def draw_disk_3d(ax, center, radius, normal=(0,0,1), color="#4f46e5", alpha=0.5, n=24):
    """Draw a disk (circle) at a point."""
    theta = np.linspace(0, 2*np.pi, n)
    normal = np.array(normal)/np.linalg.norm(normal)
    if abs(normal[2]) < 0.99: perp1 = np.cross(normal, [0,0,1])
    else: perp1 = np.cross(normal, [1,0,0])
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(normal, perp1)
    verts = []
    for t in theta:
        pt = np.array(center) + radius*(np.cos(t)*perp1 + np.sin(t)*perp2)
        verts.append(pt)
    poly = Poly3DCollection([verts], alpha=alpha, facecolor=color, edgecolor=color, linewidths=0.5)
    ax.add_collection3d(poly)


def draw_pedestal_3d(ax):
    """Draw a realistic 3-tier robot pedestal base in 3D."""
    
    draw_cylinder_3d(ax, [0,0,-BASE_H], [0,0,-BASE_H+10],
                     radius=BASE_R, color=TH["ped_side"], alpha=0.90)
    draw_disk_3d(ax, [0,0,-BASE_H], BASE_R,    color=TH["ped_top"], alpha=0.85)
    draw_disk_3d(ax, [0,0,-BASE_H+10], BASE_R, color=TH["ped_top"], alpha=0.7)
    
    t = np.linspace(0, 2*np.pi, 48)
    z0 = -BASE_H + 10
    ax.plot(BASE_R*np.cos(t), BASE_R*np.sin(t), z0,
            color=TH["ped_edge"], lw=1.2, alpha=0.7)

    
    draw_cylinder_3d(ax, [0,0,-BASE_H+10], [0,0,-PLATE_H],
                     radius=NECK_R+12, color=TH["ped_side"], alpha=0.88)
    draw_disk_3d(ax, [0,0,-PLATE_H], NECK_R+12, color=TH["ped_top"], alpha=0.65)
    
    z1, z2 = -BASE_H+10, -PLATE_H
    for angle in np.linspace(0, 2*np.pi, 7)[:-1]:
        xp = (NECK_R+12)*np.cos(angle); yp = (NECK_R+12)*np.sin(angle)
        ax.plot([xp,xp],[yp,yp],[z1,z2], color=TH["ped_edge"], lw=0.6, alpha=0.35)

    
    draw_cylinder_3d(ax, [0,0,-PLATE_H], [0,0,0],
                     radius=NECK_R, color=TH["mount_c"], alpha=0.92)
    draw_disk_3d(ax, [0,0,0],       NECK_R,   color=TH["mount_hi"], alpha=0.80)
    draw_disk_3d(ax, [0,0,-PLATE_H], NECK_R,  color=TH["ped_top"],  alpha=0.60)

    
    ring_t = np.linspace(0, 2*np.pi, 80)
    ax.plot(NECK_R*np.cos(ring_t), NECK_R*np.sin(ring_t), 0,
            color=TH["ped_ring"], lw=1.8, alpha=0.9)
    ax.plot((NECK_R+6)*np.cos(ring_t), (NECK_R+6)*np.sin(ring_t), 0,
            color=TH["ped_ring"], lw=0.5, alpha=0.3)

    
    for angle in [45, 135, 225, 315]:
        bx = (NECK_R-8)*np.cos(np.radians(angle))
        by = (NECK_R-8)*np.sin(np.radians(angle))
        draw_sphere_3d(ax, [bx, by, 1], radius=2.8, color="#94a3b8", alpha=0.9)

    
    ct = np.linspace(0, 2*np.pi, 20)
    cx, cy = (NECK_R+14)*np.cos(np.radians(270)), (NECK_R+14)*np.sin(np.radians(270))
    ax.plot([0, cx], [0, cy], [-BASE_H*0.5, -BASE_H*0.5],
            color=TH["ped_edge"], lw=2, alpha=0.25)


def draw_pedestal_side(ax, x_offset=0):
    """Draw pedestal silhouette for Side (RZ) view."""
    
    ax.fill_betweenx([-BASE_H, -BASE_H+10],
                     x_offset - BASE_R, x_offset + BASE_R,
                     color=TH["ped_side"], alpha=0.55)
    ax.plot([x_offset-BASE_R, x_offset-BASE_R], [-BASE_H, -BASE_H+10],
            color=TH["ped_edge"], lw=1.2, alpha=0.7)
    ax.plot([x_offset+BASE_R, x_offset+BASE_R], [-BASE_H, -BASE_H+10],
            color=TH["ped_edge"], lw=1.2, alpha=0.7)
    
    nw = NECK_R + 12
    ax.fill_betweenx([-BASE_H+10, -PLATE_H],
                     x_offset - nw, x_offset + nw,
                     color=TH["ped_side"], alpha=0.50)
    ax.plot([x_offset-nw, x_offset-nw], [-BASE_H+10, -PLATE_H],
            color=TH["ped_edge"], lw=1.0, alpha=0.6)
    ax.plot([x_offset+nw, x_offset+nw], [-BASE_H+10, -PLATE_H],
            color=TH["ped_edge"], lw=1.0, alpha=0.6)
    
    ax.fill_betweenx([-PLATE_H, 0],
                     x_offset - NECK_R, x_offset + NECK_R,
                     color=TH["mount_c"], alpha=0.60)
    ax.plot([x_offset-NECK_R, x_offset-NECK_R], [-PLATE_H, 0],
            color=TH["mount_hi"], lw=1.2, alpha=0.55)
    ax.plot([x_offset+NECK_R, x_offset+NECK_R], [-PLATE_H, 0],
            color=TH["mount_hi"], lw=1.2, alpha=0.55)
    # Top glow line
    ax.plot([x_offset-NECK_R, x_offset+NECK_R], [0, 0],
            color=TH["ped_ring"], lw=2.0, alpha=0.85)
    # Base shadow
    ax.plot([x_offset-BASE_R, x_offset+BASE_R], [-BASE_H, -BASE_H],
            color="#000000", lw=3, alpha=0.4)


def draw_pedestal_top(ax):
    """Draw pedestal rings for Top (XY) view."""
    tw = np.linspace(0, 2*np.pi, 80)
    
    ax.fill_between(BASE_R*np.cos(tw), BASE_R*np.sin(tw),
                    alpha=0.12, color=TH["ped_side"])
    ax.plot(BASE_R*np.cos(tw), BASE_R*np.sin(tw),
            color=TH["ped_edge"], lw=1.5, alpha=0.55, ls='--')
    
    ax.fill_between(NECK_R*np.cos(tw), NECK_R*np.sin(tw),
                    alpha=0.25, color=TH["mount_c"])
    ax.plot(NECK_R*np.cos(tw), NECK_R*np.sin(tw),
            color=TH["ped_ring"], lw=2.0, alpha=0.85)
    
    for angle in [45, 135, 225, 315]:
        bx = (NECK_R-8)*np.cos(np.radians(angle))
        by = (NECK_R-8)*np.sin(np.radians(angle))
        ax.scatter(bx, by, color="#94a3b8", s=18, zorder=6, edgecolors=TH["ped_edge"], linewidths=0.8)



class App:
    def __init__(self, root):
        self.root = root
        self.root.title("SCARA 4-DOF — Động Học Thuận & Nghịch  |  v4.0")
        self.root.configure(bg=TH["bg"])
        self.root.state("zoomed")
        self.root.minsize(1300,720)
        self._lock=False

        self.F = {
            "title":  tkfont.Font(family="Segoe UI",size=13,weight="bold"),
            "h2":     tkfont.Font(family="Segoe UI",size=10,weight="bold"),
            "body":   tkfont.Font(family="Segoe UI",size=9),
            "mono":   tkfont.Font(family="Consolas",size=9),
            "monob":  tkfont.Font(family="Consolas",size=10,weight="bold"),
            "monos":  tkfont.Font(family="Consolas",size=8),
            "badge":  tkfont.Font(family="Segoe UI",size=8,weight="bold"),
            "btn":    tkfont.Font(family="Consolas",size=10,weight="bold"),
            "mat":    tkfont.Font(family="Consolas",size=9),
            "mlab":   tkfont.Font(family="Segoe UI",size=7,weight="bold"),
            "tiny":   tkfont.Font(family="Consolas",size=7),
        }
        
        self.t1=tk.DoubleVar(value=30); self.t2=tk.DoubleVar(value=-45)
        self.d3=tk.DoubleVar(value=50); self.t4=tk.DoubleVar(value=0)
        self.l1=tk.DoubleVar(value=250); self.l2=tk.DoubleVar(value=200)
        self.ik_x=tk.DoubleVar(value=0); self.ik_y=tk.DoubleVar(value=0)
        self.ik_z=tk.DoubleVar(value=0); self.ik_phi=tk.DoubleVar(value=0)
        self.ik_cfg=tk.StringVar(value="elbow_up"); self.ik_sols={}
        self.view=tk.StringVar(value="3d")
        self.grid_on=tk.BooleanVar(value=True)
        self.trace_on=tk.BooleanVar(value=False)
        self.trace_pts=[]; self.animating=False

        
        self.traj_type=tk.StringVar(value="line")
        self.traj_lx1=tk.DoubleVar(value=300); self.traj_ly1=tk.DoubleVar(value=0)
        self.traj_lx2=tk.DoubleVar(value=200); self.traj_ly2=tk.DoubleVar(value=250)
        self.traj_lz=tk.DoubleVar(value=0)
        self.traj_cx=tk.DoubleVar(value=250); self.traj_cy=tk.DoubleVar(value=0)
        self.traj_cr=tk.DoubleVar(value=100); self.traj_cz=tk.DoubleVar(value=0)
        self.traj_N=tk.IntVar(value=150); self.traj_T=tk.DoubleVar(value=4.0)
        self.traj_phi=tk.DoubleVar(value=0)
        self.traj_running=False; self.traj_pts_planned=None
        self.traj_joints_log=[]  

        self._build()
        for v in [self.t1,self.t2,self.d3,self.t4,self.l1,self.l2]:
            v.trace_add("write", lambda *_: self._sched())
        self._sched()

    def _sched(self):
        if not self._lock:
            self._lock=True
            self.root.after(30, self._do_update)

    def _do_update(self):
        self._lock=False
        try:
            fk = forward_kinematics(self.t1.get(),self.t2.get(),self.d3.get(),
                                     self.t4.get(),self.l1.get(),self.l2.get())
        except: return
        for k,v in [("X",fk['x']),("Y",fk['y']),("Z",fk['z']),("φ",fk['phi'])]:
            self.res_l[k].config(text=f"{v:.2f}")
        t1,t2,t4=self.t1.get(),self.t2.get(),self.t4.get()
        d3v,l1v,l2v=self.d3.get(),self.l1.get(),self.l2.get()
        dh=[[f"{t1:.1f}","0",f"{l1v:.0f}","0"],
            [f"{t2:.1f}","180",f"{l2v:.0f}","0"],
            ["0","0","0",f"{d3v:.1f}"],
            [f"{t4:.1f}","0","0","0"]]
        for i in range(4):
            for j in range(4):
                self.dh_c[(i,j+1)].config(text=dh[i][j])
        T=fk["T"]
        for i in range(4):
            for j in range(4):
                v=T[i,j]
                self.mat_c[(i,j)].config(text=f"{v: .4f}",
                    fg=TH["cyan"] if abs(v)>0.001 else TH["text3"])
        wm=l1v+l2v; wn=abs(l1v-l2v); cur=np.sqrt(fk['x']**2+fk['y']**2)
        self.ws_l["max"].config(text=f"{wm:.0f}")
        self.ws_l["min"].config(text=f"{wn:.0f}")
        self.ws_l["cur"].config(text=f"{cur:.1f}")
        
        if hasattr(self, '_sb_items'):
            self._sb_items["tcp_xyz"].config(
                text=f"TCP: ({fk['x']:.0f}, {fk['y']:.0f}, {fk['z']:.0f}) mm")
            self._sb_items["reach"].config(
                text=f"Reach: {cur:.0f} / {wm:.0f} mm  [{cur/wm*100:.0f}%]")
        if self.trace_on.get():
            self.trace_pts.append((fk['x'],fk['y'],fk['z']))
            if len(self.trace_pts)>2000: self.trace_pts=self.trace_pts[-1500:]
        self._render(fk)

    
    def _build(self):
        
        hd = tk.Frame(self.root, bg=TH["bg2"], height=56)
        hd.pack(fill="x"); hd.pack_propagate(False)

        
        tk.Frame(hd, bg=TH["hi"], width=4).pack(side="left", fill="y")

        lf = tk.Frame(hd, bg=TH["bg2"]); lf.pack(side="left", padx=14, fill="y")
        
        ic = tk.Frame(lf, bg="#1e1b4b", highlightthickness=1,
                      highlightbackground=TH["purple"])
        ic.pack(side="left", padx=(0,12), pady=8)
        tk.Label(ic, text=" ⬡ ", font=("Segoe UI", 16),
                 bg="#1e1b4b", fg=TH["purple"]).pack(padx=4, pady=2)

        ttf = tk.Frame(lf, bg=TH["bg2"]); ttf.pack(side="left", fill="y", pady=6)
        tk.Label(ttf, text="SCARA Robot Controller",
                 font=tkfont.Font(family="Segoe UI", size=14, weight="bold"),
                 bg=TH["bg2"], fg=TH["text"]).pack(anchor="w")
        sub_f = tk.Frame(ttf, bg=TH["bg2"]); sub_f.pack(anchor="w")
        tk.Label(sub_f, text="Forward & Inverse Kinematics",
                 font=self.F["body"], bg=TH["bg2"], fg=TH["text3"]).pack(side="left")
        tk.Label(sub_f, text=" · ", font=self.F["body"],
                 bg=TH["bg2"], fg=TH["border"]).pack(side="left")
        tk.Label(sub_f, text="4 DOF",
                 font=tkfont.Font(family="Consolas", size=9, weight="bold"),
                 bg=TH["bg2"], fg=TH["cyan"]).pack(side="left")
        tk.Label(sub_f, text=" · Base + 3 Links + Gripper",
                 font=self.F["body"], bg=TH["bg2"], fg=TH["text3"]).pack(side="left")

        rf = tk.Frame(hd, bg=TH["bg2"]); rf.pack(side="right", padx=14, fill="y")

        def _mk_hbtn(parent, text, fg, cmd):
            f = tk.Frame(parent, bg=TH["card2"], highlightthickness=1,
                         highlightbackground=TH["border"])
            f.pack(side="right", padx=4, pady=10)
            b = tk.Button(f, text=text, font=self.F["badge"],
                bg=TH["card2"], fg=fg, bd=0, padx=16, pady=5,
                cursor="hand2", command=cmd, activebackground=TH["input"],
                activeforeground=fg)
            b.pack()
            return b

        self.btn_an = _mk_hbtn(rf, "▶  Demo", TH["cyan"], self._toggle_anim)
        _mk_hbtn(rf, "↺  Reset", TH["rose"], self._reset)

        main = tk.Frame(self.root, bg=TH["bg"])
        main.pack(fill="both", expand=True, padx=8, pady=(4,0))
        main.columnconfigure(0, weight=0, minsize=330)
        main.columnconfigure(1, weight=1)
        main.columnconfigure(2, weight=0, minsize=350)
        main.rowconfigure(0, weight=1)
        self._build_left(main)
        self._build_center(main)
        self._build_right(main)

        
        sb = tk.Frame(self.root, bg=TH["bg2"], height=24)
        sb.pack(fill="x"); sb.pack_propagate(False)
        tk.Frame(sb, bg=TH["hi"], width=4).pack(side="left", fill="y")
        self._sb_items = {}
        for key, label, col in [
            ("tcp_xyz", "TCP: —", TH["text3"]),
            ("base_h",  f"Base H: {BASE_H} mm", TH["teal"]),
            ("reach",   "Reach: —", TH["cyan"]),
            ("mode",    "Mode: FK", TH["purple"]),
        ]:
            lbl = tk.Label(sb, text=label, font=self.F["tiny"],
                           bg=TH["bg2"], fg=col)
            lbl.pack(side="left", padx=12)
            self._sb_items[key] = lbl
        tk.Label(sb, text="SCARA 4-DOF v4.0", font=self.F["tiny"],
                 bg=TH["bg2"], fg=TH["text3"]).pack(side="right", padx=12)

    def _build_left(self, par):
        o = tk.Frame(par, bg=TH["card"], highlightthickness=1,
                     highlightbackground=TH["border"])
        o.grid(row=0, column=0, sticky="nsew", padx=(0,4))
        cv = tk.Canvas(o, bg=TH["card"], highlightthickness=0, bd=0)
        sb = tk.Scrollbar(o, orient="vertical", command=cv.yview)
        p = tk.Frame(cv, bg=TH["card"], padx=12, pady=8)
        p.bind("<Configure>", lambda e: cv.configure(scrollregion=cv.bbox("all")))
        cv.create_window((0,0), window=p, anchor="nw")
        cv.configure(yscrollcommand=sb.set)
        cv.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        cv.bind_all("<MouseWheel>", lambda e: cv.yview_scroll(int(-1*(e.delta/120)),"units"), add="+")

        
        self._hdr(p, "Điều Khiển Khớp", "FK", TH["purple"])
        for nm,var,lo,hi,u,col in [
            ("θ₁  Khớp 1",self.t1,-150,150,"°",TH["j1"]),
            ("θ₂  Khớp 2",self.t2,-150,150,"°",TH["j2"]),
            ("d₃  Khớp 3",self.d3,0,150,"mm",TH["j3"]),
            ("θ₄  Khớp 4",self.t4,-360,360,"°",TH["j4"])]:
            self._slider(p, nm, var, lo, hi, u, col)

        self._sep(p)
        self._hdr(p, "Thông Số Link", "DH", TH["cyan"])
        for nm,var in [("l₁ — Link 1",self.l1),("l₂ — Link 2",self.l2)]:
            self._param_row(p, nm, var, "mm")

        
        self._sep(p, TH["pink"])
        self._hdr(p, "Động Học Nghịch", "IK", TH["pink"])
        for nm,var,u,col in [
            ("X đích",self.ik_x,"mm",TH["j1"]),("Y đích",self.ik_y,"mm",TH["j2"]),
            ("Z đích",self.ik_z,"mm",TH["j3"]),("φ đích",self.ik_phi,"°",TH["j4"])]:
            self._ik_row(p, nm, var, u, col)

        cf = tk.Frame(p, bg=TH["card"]); cf.pack(fill="x", pady=(6,4))
        tk.Label(cf, text="Cấu hình:", font=self.F["body"], bg=TH["card"],
                 fg=TH["text2"]).pack(side="left")
        for txt,val,col in [("⬆ Up","elbow_up",TH["green"]),
                             ("⬇ Down","elbow_down",TH["purple"])]:
            tk.Radiobutton(cf, text=txt, variable=self.ik_cfg, value=val,
                           font=self.F["badge"], bg=TH["card"], fg=col,
                           selectcolor=TH["input"], activebackground=TH["card"],
                           activeforeground=col).pack(side="left", padx=5)

        bf = tk.Frame(p, bg=TH["card"]); bf.pack(fill="x", pady=(6,3))
        tk.Button(bf, text="📐 Tính IK", font=self.F["badge"],
                  bg="#7c3aed", fg="white", bd=0, padx=14, pady=7,
                  cursor="hand2", command=self._calc_ik
                  ).pack(side="left", padx=(0,3), expand=True, fill="x")
        tk.Button(bf, text="✅ Áp dụng", font=self.F["badge"],
                  bg="#059669", fg="white", bd=0, padx=14, pady=7,
                  cursor="hand2", command=self._apply_ik
                  ).pack(side="left", padx=(3,0), expand=True, fill="x")
        tk.Button(p, text="📍 Lấy vị trí hiện tại → IK", font=self.F["badge"],
                  bg=TH["input"], fg=TH["cyan"], bd=0, padx=10, pady=5,
                  cursor="hand2", command=self._fill_ik).pack(fill="x", pady=3)

        self.ik_st = tk.Label(p, text="", font=self.F["body"], bg=TH["card"],
                              fg=TH["text3"], wraplength=280, justify="left")
        self.ik_st.pack(fill="x")
        self.ik_sf = tk.Frame(p, bg=TH["card"]); self.ik_sf.pack(fill="x")

        
        self._sep(p, TH["cyan"])
        self._hdr(p, "Lập Quỹ Đạo", "TRAJ", TH["cyan"])

        
        tyf = tk.Frame(p, bg=TH["card"]); tyf.pack(fill="x", pady=(0,6))
        tk.Label(tyf, text="Loại:", font=self.F["body"], bg=TH["card"],
                 fg=TH["text2"]).pack(side="left")
        for txt,val,col in [("📏 Thẳng","line",TH["green"]),
                             ("⭕ Tròn","circle",TH["cyan"])]:
            tk.Radiobutton(tyf, text=txt, variable=self.traj_type, value=val,
                           font=self.F["badge"], bg=TH["card"], fg=col,
                           selectcolor=TH["input"], activebackground=TH["card"],
                           activeforeground=col,
                           command=self._traj_type_changed).pack(side="left", padx=5)

        
        self.traj_line_f = tk.Frame(p, bg=TH["card"])
        self.traj_line_f.pack(fill="x")
        for nm,var,u in [("X₁ start",self.traj_lx1,"mm"),("Y₁ start",self.traj_ly1,"mm"),
                          ("X₂ end",self.traj_lx2,"mm"),("Y₂ end",self.traj_ly2,"mm"),
                          ("Z",self.traj_lz,"mm")]:
            self._ik_row(self.traj_line_f, nm, var, u, TH["green"])

        
        self.traj_circ_f = tk.Frame(p, bg=TH["card"])
        for nm,var,u in [("Tâm X",self.traj_cx,"mm"),("Tâm Y",self.traj_cy,"mm"),
                          ("Bán kính",self.traj_cr,"mm"),("Z",self.traj_cz,"mm")]:
            self._ik_row(self.traj_circ_f, nm, var, u, TH["cyan"])

        
        comf = tk.Frame(p, bg=TH["card"]); comf.pack(fill="x", pady=(4,0))
        for nm,var,u in [("Số điểm N",self.traj_N,""),("Thời gian T",self.traj_T,"s"),
                          ("φ hướng",self.traj_phi,"°")]:
            self._ik_row(comf, nm, var, u, TH["text2"])

        
        tbf = tk.Frame(p, bg=TH["card"]); tbf.pack(fill="x", pady=(8,3))
        self.btn_traj = tk.Button(tbf, text="▶  Chạy quỹ đạo", font=self.F["badge"],
                  bg="#059669", fg="white", bd=0, padx=14, pady=7,
                  cursor="hand2", command=self._run_trajectory)
        self.btn_traj.pack(side="left", padx=(0,3), expand=True, fill="x")
        tk.Button(tbf, text="📊 Đồ thị", font=self.F["badge"],
                  bg=TH["input"], fg=TH["cyan"], bd=0, padx=14, pady=7,
                  cursor="hand2", command=self._show_traj_plots
                  ).pack(side="left", padx=(3,0), expand=True, fill="x")

        self.traj_st = tk.Label(p, text="", font=self.F["body"], bg=TH["card"],
                                fg=TH["text3"], wraplength=280, justify="left")
        self.traj_st.pack(fill="x")

    def _build_center(self, par):
        panel = tk.Frame(par, bg=TH["card"], highlightthickness=1,
                         highlightbackground=TH["border"])
        panel.grid(row=0, column=1, sticky="nsew", padx=4)

        tf = tk.Frame(panel, bg=TH["card2"], padx=4, pady=4)
        tf.pack(fill="x", pady=(4,2), padx=4)
        self.tabs = {}
        for v,l,ico in [("3d","3D","🔮"),("top","Top XY","🔵"),("side","Side RZ","🟩")]:
            b = tk.Button(tf, text=f"{ico} {l}", font=self.F["badge"],
                bg=TH["hi"] if v=="3d" else TH["card"],
                fg=TH["text"] if v=="3d" else TH["text3"],
                bd=0, padx=16, pady=6, cursor="hand2",
                command=lambda v=v: self._set_view(v))
            b.pack(side="left", padx=2, expand=True, fill="x")
            self.tabs[v] = b

        of = tk.Frame(panel, bg=TH["card"]); of.pack(fill="x", padx=6)
        for txt,var in [("Grid",self.grid_on),("Trace",self.trace_on)]:
            tk.Checkbutton(of, text=txt, variable=var, font=self.F["badge"],
                           bg=TH["card"], fg=TH["text2"], selectcolor=TH["input"],
                           activebackground=TH["card"],
                           command=self._sched).pack(side="left", padx=6)
        tk.Button(of, text="Xóa Trace", font=self.F["badge"], bg=TH["card"],
                  fg=TH["text3"], bd=0, padx=8, cursor="hand2",
                  command=self._clear_trace).pack(side="left", padx=6)

        self.fig = Figure(facecolor=TH["card"], dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=panel)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=4, pady=4)

    def _build_right(self, par):
        o = tk.Frame(par, bg=TH["card"], highlightthickness=1,
                     highlightbackground=TH["border"])
        o.grid(row=0, column=2, sticky="nsew", padx=(4,0))
        cv = tk.Canvas(o, bg=TH["card"], highlightthickness=0, bd=0)
        sb = tk.Scrollbar(o, orient="vertical", command=cv.yview)
        p = tk.Frame(cv, bg=TH["card"], padx=12, pady=8)
        p.bind("<Configure>", lambda e: cv.configure(scrollregion=cv.bbox("all")))
        cv.create_window((0,0), window=p, anchor="nw")
        cv.configure(yscrollcommand=sb.set)
        cv.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

        
        self._hdr(p, "Vị Trí TCP", "XYZφ", TH["rose"])
        self.res_l = {}
        rg = tk.Frame(p, bg=TH["card"]); rg.pack(fill="x", pady=(0,6))
        rg.columnconfigure(0, weight=1); rg.columnconfigure(1, weight=1)
        for i,(nm,col) in enumerate([("X",TH["j1"]),("Y",TH["j2"]),
                                      ("Z",TH["j3"]),("φ",TH["j4"])]):
            cd = tk.Frame(rg, bg=TH["input"], highlightthickness=1,
                          highlightbackground=TH["border"], padx=8, pady=5)
            cd.grid(row=i//2, column=i%2, padx=2, pady=2, sticky="nsew")
            tk.Label(cd, text=nm, font=self.F["monob"], bg=TH["input"],
                     fg=col).pack(side="left")
            u = "°" if nm=="φ" else "mm"
            vf = tk.Frame(cd, bg=TH["input"]); vf.pack(side="right")
            lbl = tk.Label(vf, text="0.00", font=self.F["monob"],
                           bg=TH["input"], fg=TH["text"])
            lbl.pack(anchor="e")
            tk.Label(vf, text=u, font=self.F["tiny"], bg=TH["input"],
                     fg=TH["text3"]).pack(anchor="e")
            self.res_l[nm] = lbl

        
        self._sep(p)
        self._hdr(p, "Bảng DH Parameters", "D-H", TH["cyan"])
        tf = tk.Frame(p, bg=TH["input"], highlightthickness=1,
                      highlightbackground=TH["border"])
        tf.pack(fill="x", pady=(0,6))
        for j,h in enumerate(["Khớp","θᵢ","αᵢ","aᵢ","dᵢ"]):
            tk.Label(tf, text=h, font=self.F["badge"], bg="#181b48",
                     fg=TH["purple"], padx=5, pady=4).grid(row=0,column=j,sticky="nsew")
            tf.columnconfigure(j, weight=1)
        self.dh_c={}
        jcol=[TH["j1"],TH["j2"],TH["j3"],TH["j4"]]
        for i in range(4):
            for j in range(5):
                c=tk.Label(tf, text="—", font=self.F["monos"], bg=TH["input"],
                           fg=TH["text2"], padx=3, pady=3)
                c.grid(row=i+1, column=j, sticky="nsew")
                self.dh_c[(i,j)]=c
            self.dh_c[(i,0)].config(text=f" {i+1}", fg=jcol[i], font=self.F["badge"])

        
        self._sep(p)
        self._hdr(p, "Ma Trận T₀₄", "n o a p", TH["cyan"])
        mf = tk.Frame(p, bg=TH["card"])
        mf.pack(fill="x", pady=(0,6))
        self.mat_c = {}
        self._matrix_widget(mf, self.mat_c, TH["cyan"])

        
        self._sep(p)
        self._hdr(p, "Không Gian Làm Việc", "WS", TH["teal"])
        self.ws_l = {}
        for k,t in [("max","R max"),("min","R min"),("cur","R hiện tại")]:
            r = tk.Frame(p, bg=TH["input"], padx=8, pady=4)
            r.pack(fill="x", pady=1)
            tk.Label(r, text=t, font=self.F["body"], bg=TH["input"],
                     fg=TH["text3"]).pack(side="left")
            v = tk.Label(r, text="—", font=self.F["mono"], bg=TH["input"],
                         fg=TH["text"])
            v.pack(side="right")
            self.ws_l[k] = v

        
        self._sep(p)
        self._hdr(p, "Thông Số Đế Robot", "BASE", TH["ped_ring"])
        base_info = [
            ("Chiều cao đế",    f"{BASE_H} mm",      TH["text"]),
            ("Bán kính chân",   f"{BASE_R} mm",       TH["text"]),
            ("Bán kính cổ",     f"{NECK_R} mm",       TH["text"]),
            ("Tấm lắp đặt",     f"{PLATE_H} mm",      TH["text"]),
            ("Số bu-lông",       "4 × M8",             TH["amber"]),
        ]
        for label, val, col in base_info:
            r = tk.Frame(p, bg=TH["input"], padx=8, pady=3,
                         highlightthickness=1, highlightbackground=TH["border"])
            r.pack(fill="x", pady=1)
            tk.Label(r, text=label, font=self.F["body"],
                     bg=TH["input"], fg=TH["text3"]).pack(side="left")
            tk.Label(r, text=val, font=self.F["mono"],
                     bg=TH["input"], fg=col).pack(side="right")

    
    def _hdr(self, p, title, badge, col):
        f=tk.Frame(p, bg=TH["card"]); f.pack(fill="x", pady=(2,6))
        # Accent bar
        tk.Frame(f, bg=col, width=3, height=16).pack(side="left", padx=(0,8))
        tk.Label(f, text=title, font=self.F["h2"], bg=TH["card"],
                 fg=TH["text"]).pack(side="left")
        tk.Label(f, text=badge, font=self.F["badge"], bg=TH["card"],
                 fg=col).pack(side="right")

    def _sep(self, p, col=None):
        c = col or TH["border"]
        tk.Frame(p, bg=c, height=1 if not col else 2).pack(fill="x", pady=12)

    def _slider(self, p, nm, var, lo, hi, u, col):
        cd = tk.Frame(p, bg=TH["input"], highlightthickness=1,
                      highlightbackground=TH["border"], padx=10, pady=7)
        cd.pack(fill="x", pady=3)
        top = tk.Frame(cd, bg=TH["input"]); top.pack(fill="x")
        d = tk.Canvas(top, width=8, height=8, bg=TH["input"], highlightthickness=0)
        d.pack(side="left", padx=(0,5)); d.create_oval(0,0,8,8, fill=col, outline=col)
        tk.Label(top, text=nm, font=self.F["body"], bg=TH["input"],
                 fg=TH["text"]).pack(side="left")
        vf = tk.Frame(top, bg="#090b18", highlightthickness=1,
                      highlightbackground=TH["border"], padx=5, pady=1)
        vf.pack(side="right")
        tk.Entry(vf, textvariable=var, font=self.F["mono"], bg="#090b18",
                 fg=col, insertbackground=TH["text"], bd=0, width=7,
                 justify="right").pack(side="left")
        tk.Label(vf, text=u, font=self.F["tiny"], bg="#090b18",
                 fg=TH["text3"]).pack(side="left", padx=(2,0))
        mid = tk.Frame(cd, bg=TH["input"]); mid.pack(fill="x", pady=(5,1))
        tk.Button(mid, text="−", font=self.F["btn"], bg=TH["trough"], fg=col,
                  activebackground=col, activeforeground="white", bd=0, width=3,
                  repeatdelay=300, repeatinterval=50, cursor="hand2",
                  command=lambda: self._adj(var,-1,lo,hi)).pack(side="left")
        tk.Scale(mid, from_=lo, to=hi, orient="horizontal", variable=var,
                 resolution=0.5, showvalue=False, bg=TH["input"], fg=col,
                 troughcolor=TH["trough"], activebackground=col,
                 highlightthickness=0, sliderrelief="flat", sliderlength=16,
                 width=12, bd=0).pack(side="left", fill="x", expand=True, padx=4)
        tk.Button(mid, text="+", font=self.F["btn"], bg=TH["trough"], fg=col,
                  activebackground=col, activeforeground="white", bd=0, width=3,
                  repeatdelay=300, repeatinterval=50, cursor="hand2",
                  command=lambda: self._adj(var,1,lo,hi)).pack(side="right")

    def _adj(self, var, d, lo, hi):
        try: var.set(max(lo,min(hi,var.get()+d)))
        except: pass

    def _param_row(self, p, nm, var, u):
        r=tk.Frame(p, bg=TH["card"]); r.pack(fill="x", pady=2)
        tk.Label(r, text=nm, font=self.F["body"], bg=TH["card"],
                 fg=TH["text2"]).pack(side="left")
        vf=tk.Frame(r, bg="#090b18", highlightthickness=1,
                    highlightbackground=TH["border"], padx=5, pady=1)
        vf.pack(side="right")
        tk.Entry(vf, textvariable=var, font=self.F["mono"], bg="#090b18",
                 fg=TH["text"], insertbackground=TH["text"], bd=0, width=6,
                 justify="right").pack(side="left")
        tk.Label(vf, text=u, font=self.F["tiny"], bg="#090b18",
                 fg=TH["text3"]).pack(side="left", padx=(2,0))

    def _ik_row(self, p, nm, var, u, col):
        r=tk.Frame(p, bg=TH["input"], highlightthickness=1,
                   highlightbackground=TH["border"], padx=8, pady=4)
        r.pack(fill="x", pady=2)
        d=tk.Canvas(r, width=8, height=8, bg=TH["input"], highlightthickness=0)
        d.pack(side="left", padx=(0,5)); d.create_oval(0,0,8,8, fill=col, outline=col)
        tk.Label(r, text=nm, font=self.F["body"], bg=TH["input"],
                 fg=TH["text"]).pack(side="left")
        vf=tk.Frame(r, bg="#090b18", highlightthickness=1,
                    highlightbackground=TH["border"], padx=5, pady=1)
        vf.pack(side="right")
        tk.Entry(vf, textvariable=var, font=self.F["mono"], bg="#090b18",
                 fg=col, insertbackground=TH["text"], bd=0, width=8,
                 justify="right").pack(side="left")
        tk.Label(vf, text=u, font=self.F["tiny"], bg="#090b18",
                 fg=TH["text3"]).pack(side="left", padx=(2,0))

    def _matrix_widget(self, parent, cells, col):
        """4×4 matrix with n,o,a,p headers and x,y,z row labels."""
        mf=tk.Frame(parent, bg="#080a14", highlightthickness=1,
                    highlightbackground=TH["border"], padx=8, pady=6)
        mf.pack(fill="x")
       
        eq = tk.Frame(mf, bg="#080a14"); eq.pack(fill="x", pady=(0,4))
        tk.Label(eq, text="T₀₄ =", font=self.F["monob"], bg="#080a14",
                 fg=col).pack(side="left")
        
        hf=tk.Frame(mf, bg="#080a14"); hf.pack(fill="x")
        tk.Label(hf, text="", width=3, bg="#080a14").pack(side="left")
        for h in ["nᵢ","oᵢ","aᵢ","pᵢ"]:
            tk.Label(hf, text=h, font=self.F["mlab"], bg="#080a14",
                     fg=col, width=9, anchor="center").pack(side="left", expand=True)
        rl = ["x","y","z",""]
        for i in range(4):
            rf=tk.Frame(mf, bg="#080a14"); rf.pack(fill="x")
            lbl_fg = col if rl[i] else TH["text3"]
            tk.Label(rf, text=rl[i], font=self.F["mlab"], bg="#080a14",
                     fg=lbl_fg, width=3, anchor="center").pack(side="left")
            for j in range(4):
                cell=tk.Label(rf, text="0.0000", font=self.F["mat"], bg="#080a14",
                              fg=TH["text2"], width=9, anchor="center", pady=1)
                cell.pack(side="left", expand=True, padx=1)
                cells[(i,j)]=cell

    def _set_view(self, v):
        self.view.set(v)
        for k,b in self.tabs.items():
            b.config(bg=TH["hi"] if k==v else TH["card"],
                     fg=TH["text"] if k==v else TH["text3"])
        self._sched()

    
    def _render(self, fk):
        self.fig.clf()
        v = self.view.get()
        if v=="3d": self._r3d(fk)
        elif v=="top": self._rtop(fk)
        else: self._rside(fk)
        self.canvas.draw_idle()

    
    def _r3d(self, fk):
        ax = self.fig.add_subplot(111, projection='3d', facecolor=TH["card"])
        pos = fk["pos"]
        l1v,l2v = self.l1.get(),self.l2.get()
        wm = l1v+l2v; lim = wm*1.15

        ax.set_xlabel("X", color=TH["text3"], fontsize=7, labelpad=4)
        ax.set_ylabel("Y", color=TH["text3"], fontsize=7, labelpad=4)
        ax.set_zlabel("Z", color=TH["text3"], fontsize=7, labelpad=4)
        ax.tick_params(colors=TH["text3"], labelsize=5, pad=1)
        for pa in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
            pa.fill=False; pa.set_edgecolor(TH["grid_c"])
        ax.grid(self.grid_on.get(), color=TH["grid_c"], alpha=0.3)
        ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.set_zlim(-BASE_H-20, 60)

        
        fg_range = np.linspace(-lim, lim, 12)
        for val in fg_range:
            ax.plot([val,val],[-lim,lim],[-BASE_H,-BASE_H], color=TH["grid_c"], lw=0.3, alpha=0.3)
            ax.plot([-lim,lim],[val,val],[-BASE_H,-BASE_H], color=TH["grid_c"], lw=0.3, alpha=0.3)

        
        tw=np.linspace(0,2*np.pi,60)
        ax.plot(wm*np.cos(tw), wm*np.sin(tw), 0, color=TH["purple"], alpha=0.12, lw=1)

       
        draw_pedestal_3d(ax)

        
        draw_disk_3d(ax, [0,0,0], 35, color=TH["base_c"], alpha=0.6)
        draw_disk_3d(ax, [0,0,2], 25, color="#8183f7", alpha=0.4)

        
        draw_cylinder_3d(ax, pos[0], pos[1], radius=14, color=TH["link1_a"], alpha=0.85)

        
        draw_sphere_3d(ax, pos[1], radius=16, color=TH["j1"], alpha=0.85)

        
        draw_cylinder_3d(ax, pos[1], pos[2], radius=11, color=TH["link2_a"], alpha=0.85)

        
        draw_sphere_3d(ax, pos[2], radius=14, color=TH["j2"], alpha=0.85)

        
        if abs(pos[2][2]-pos[3][2]) > 1:
            ax.plot([pos[2][0],pos[3][0]],[pos[2][1],pos[3][1]],[pos[2][2],pos[3][2]],
                    color=TH["j3"], lw=4, ls='--', alpha=0.8, solid_capstyle='round')

        
        draw_sphere_3d(ax, pos[3], radius=12, color=TH["j3"], alpha=0.85)

        
        g = fk["grip"]; tcp = pos[-1]

        
        draw_cylinder_3d(ax, g["body_top"], g["body_bot"],
                         radius=5, color=TH["grip_c"], alpha=0.9)

        
        ax.plot([g['lb'][0],g['rb'][0]],[g['lb'][1],g['rb'][1]],[g['lb'][2],g['rb'][2]],
                color=TH["grip_c"], lw=5, solid_capstyle='round', alpha=0.95)

        
        for base, tip in [(g['lb'], g['lt']), (g['rb'], g['rt'])]:
            draw_cylinder_3d(ax, base, tip, radius=4, color=TH["grip_c"], alpha=0.92)

        
        for tip in [g['lt'], g['rt']]:
            draw_sphere_3d(ax, tip, radius=5, color=TH["amber"], alpha=0.95)

        
        draw_sphere_3d(ax, tcp, radius=9, color=TH["tcp_c"], alpha=0.9)

        
        ax.plot([tcp[0], g['ar'][0]],[tcp[1], g['ar'][1]],[tcp[2], g['ar'][2]],
                color=TH["tcp_c"], lw=2, alpha=0.8, ls='--')
        ax.scatter(*g['ar'], color=TH["tcp_c"], s=20, marker='^', zorder=10)

        
        ax.text(tcp[0], tcp[1], tcp[2]+22, f"φ={fk['phi']:.1f}°",
                color=TH["grip_c"], fontsize=7, ha='center', fontweight='bold')

        
        xs=[p[0] for p in pos]; ys=[p[1] for p in pos]
        ax.plot(xs, ys, zs=-BASE_H, color=TH["text3"], alpha=0.06, lw=2)

        
        if self.trace_on.get() and len(self.trace_pts)>1:
            tp=np.array(self.trace_pts); n=len(tp)
            for i in range(1,n):
                ax.plot(tp[i-1:i+1,0],tp[i-1:i+1,1],tp[i-1:i+1,2],
                        color=TH["tcp_c"], alpha=0.1+0.7*(i/n), lw=1)

       
        if self.traj_pts_planned is not None and len(self.traj_pts_planned)>1:
            pp=self.traj_pts_planned
            ax.plot(pp[:,0], pp[:,1], pp[:,2], color=TH["cyan"], alpha=0.3, lw=1.5, ls='--')

        
        al = 30
        ax.quiver(tcp[0],tcp[1],tcp[2], al,0,0, color="#f87171", arrow_length_ratio=0.15, lw=1.2, alpha=0.7)
        ax.quiver(tcp[0],tcp[1],tcp[2], 0,al,0, color="#4ade80", arrow_length_ratio=0.15, lw=1.2, alpha=0.7)
        ax.quiver(tcp[0],tcp[1],tcp[2], 0,0,al, color="#60a5fa", arrow_length_ratio=0.15, lw=1.2, alpha=0.7)

        ax.view_init(elev=28, azim=-55)
        ax.set_title("SCARA 4-DOF — 3D", color=TH["text2"], fontsize=9, pad=6)

    
    def _rtop(self, fk):
        ax=self.fig.add_subplot(111, facecolor=TH["card"])
        pos=fk["pos"]; l1v,l2v=self.l1.get(),self.l2.get()
        wm=l1v+l2v; wn=abs(l1v-l2v); lim=wm*1.3
        ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.set_aspect('equal')
        ax.set_xlabel("X (mm)", color=TH["text3"], fontsize=8)
        ax.set_ylabel("Y (mm)", color=TH["text3"], fontsize=8)
        ax.tick_params(colors=TH["text3"], labelsize=6)
        if self.grid_on.get(): ax.grid(True, color=TH["grid_c"], alpha=0.3, lw=0.5)
        ax.axhline(0, color=TH["grid_c"], lw=0.5, alpha=0.3)
        ax.axvline(0, color=TH["grid_c"], lw=0.5, alpha=0.3)

        tw=np.linspace(0,2*np.pi,200)
        ax.fill_between(wm*np.cos(tw),wm*np.sin(tw), alpha=0.025, color=TH["purple"])
        ax.plot(wm*np.cos(tw),wm*np.sin(tw), color=TH["purple"], alpha=0.2, lw=1, ls='--')
        if wn>5: ax.plot(wn*np.cos(tw),wn*np.sin(tw), color=TH["purple"], alpha=0.1, lw=1, ls='--')

        
        draw_pedestal_top(ax)

        xs,ys=[p[0] for p in pos],[p[1] for p in pos]

        
        ax.plot([xs[0],xs[1]],[ys[0],ys[1]], color="black", lw=14, solid_capstyle='round', alpha=0.12, zorder=1)
        ax.plot([xs[1],xs[2]],[ys[1],ys[2]], color="black", lw=11, solid_capstyle='round', alpha=0.12, zorder=1)

        
        ax.plot([xs[0],xs[1]],[ys[0],ys[1]], color=TH["link1_a"], lw=11, solid_capstyle='round', alpha=0.9, zorder=2)
        ax.plot([xs[0],xs[1]],[ys[0],ys[1]], color=TH["link1_b"], lw=5, solid_capstyle='round', alpha=0.4, zorder=3)
        ax.plot([xs[1],xs[2]],[ys[1],ys[2]], color=TH["link2_a"], lw=9, solid_capstyle='round', alpha=0.9, zorder=2)
        ax.plot([xs[1],xs[2]],[ys[1],ys[2]], color=TH["link2_b"], lw=4, solid_capstyle='round', alpha=0.4, zorder=3)

        
        jc=[TH["base_c"],TH["j1"],TH["j2"],TH["j3"],TH["tcp_c"]]
        jl=["Base","J1","J2","J3","TCP"]; jsz=[160,140,120,90,150]
        for i in [0,1,2,4]:
            ax.scatter(xs[i],ys[i], color=jc[i], s=jsz[i]*2.5, alpha=0.08, zorder=4, edgecolors='none')
            ax.scatter(xs[i],ys[i], color=jc[i], s=jsz[i]*1.3, alpha=0.15, zorder=4, edgecolors='none')
            ax.scatter(xs[i],ys[i], color=jc[i], s=jsz[i], zorder=5, edgecolors='white', linewidths=1.5)
            ax.annotate(jl[i], (xs[i],ys[i]), textcoords="offset points", xytext=(12,12),
                        fontsize=7, color=jc[i], fontweight='bold', alpha=0.9,
                        bbox=dict(boxstyle='round,pad=0.2', facecolor=TH["card"],
                                  edgecolor=jc[i], alpha=0.7))

        
        g=fk["grip"]; tcp=pos[-1]

        
        for base in [g['lb'], g['rb']]:
            ax.scatter(base[0], base[1], color=TH["grip_c"],
                       s=120, zorder=6, edgecolors='white', linewidths=1.2)
            ax.scatter(base[0], base[1], color=TH["grip_c"],
                       s=300, alpha=0.15, zorder=5, edgecolors='none')

        
        ax.plot([g['lb'][0],g['rb'][0]],[g['lb'][1],g['rb'][1]],
                color=TH["grip_c"], lw=5, solid_capstyle='round', alpha=0.9, zorder=4)

        
        ax.annotate('', xy=(g['ar'][0],g['ar'][1]), xytext=(tcp[0],tcp[1]),
                    arrowprops=dict(arrowstyle='->', color=TH["tcp_c"], lw=2.5), zorder=6)

        ax.annotate(f"φ={fk['phi']:.1f}°", xy=(g['ar'][0],g['ar'][1]),
                    textcoords='offset points', xytext=(8,5), fontsize=7,
                    color=TH["grip_c"], fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor=TH["card"],
                              edgecolor=TH["grip_c"], alpha=0.8))

        
        al=25
        ax.annotate('', xy=(tcp[0]+al,tcp[1]), xytext=(tcp[0],tcp[1]),
                    arrowprops=dict(arrowstyle='->', color='#f87171', lw=1.5), zorder=7)
        ax.annotate('', xy=(tcp[0],tcp[1]+al), xytext=(tcp[0],tcp[1]),
                    arrowprops=dict(arrowstyle='->', color='#4ade80', lw=1.5), zorder=7)
        ax.text(tcp[0]+al+3, tcp[1], "x", fontsize=6, color='#f87171')
        ax.text(tcp[0], tcp[1]+al+3, "y", fontsize=6, color='#4ade80')

        
        if self.trace_on.get() and len(self.trace_pts)>1:
            tp=np.array(self.trace_pts); n=len(tp)
            for i in range(1,n):
                ax.plot(tp[i-1:i+1,0],tp[i-1:i+1,1],
                        color=TH["tcp_c"], alpha=0.05+0.8*(i/n), lw=1.5)

        
        if self.traj_pts_planned is not None and len(self.traj_pts_planned)>1:
            pp=self.traj_pts_planned
            ax.plot(pp[:,0], pp[:,1], color=TH["cyan"], alpha=0.3, lw=1.5, ls='--', zorder=1)
            ax.plot(pp[0,0], pp[0,1], 'o', color=TH["green"], ms=6, alpha=0.6, zorder=8)
            ax.plot(pp[-1,0], pp[-1,1], 's', color=TH["rose"], ms=6, alpha=0.6, zorder=8)

        ax.annotate(f"({fk['x']:.0f}, {fk['y']:.0f})", xy=(tcp[0],tcp[1]),
                    textcoords='offset points', xytext=(15,-22), fontsize=7,
                    color=TH["tcp_c"],
                    bbox=dict(boxstyle='round,pad=0.3', facecolor=TH["card"],
                              edgecolor=TH["tcp_c"], alpha=0.85))
        ax.set_title("SCARA — Top View (XY)", color=TH["text2"], fontsize=9, pad=6)

    
    def _rside(self, fk):
        ax=self.fig.add_subplot(111, facecolor=TH["card"])
        pos=fk["pos"]; l1v,l2v=self.l1.get(),self.l2.get()
        wm=l1v+l2v
        radii=[np.sqrt(p[0]**2+p[1]**2) for p in pos]
        zs=[p[2] for p in pos]
        lx=wm*1.3
        ax.set_xlim(-lx*0.25, lx); ax.set_ylim(-BASE_H-30, 50); ax.set_aspect('equal')
        ax.set_xlabel("R (mm)", color=TH["text3"], fontsize=8)
        ax.set_ylabel("Z (mm)", color=TH["text3"], fontsize=8)
        ax.tick_params(colors=TH["text3"], labelsize=6)
        if self.grid_on.get(): ax.grid(True, color=TH["grid_c"], alpha=0.3, lw=0.5)
        ax.axhline(0, color=TH["grid_c"], lw=0.8, alpha=0.4)

        
        draw_pedestal_side(ax, x_offset=0)

        
        ax.axhline(-BASE_H, color=TH["floor"], lw=1.5, alpha=0.5, ls='-')
        ax.fill_between([-lx*0.25, lx], -BASE_H-30, -BASE_H,
                         color=TH["floor"], alpha=0.4)

        
        ax.plot([0,radii[1]],[zs[0],zs[1]], color="black", lw=14, solid_capstyle='round', alpha=0.12, zorder=1)
        ax.plot([radii[1],radii[2]],[zs[1],zs[2]], color="black", lw=11, solid_capstyle='round', alpha=0.12, zorder=1)

       
        ax.plot([0,radii[1]],[zs[0],zs[1]], color=TH["link1_a"], lw=11, solid_capstyle='round', alpha=0.9, zorder=2)
        ax.plot([0,radii[1]],[zs[0],zs[1]], color=TH["link1_b"], lw=5, solid_capstyle='round', alpha=0.4, zorder=3)
        ax.plot([radii[1],radii[2]],[zs[1],zs[2]], color=TH["link2_a"], lw=9, solid_capstyle='round', alpha=0.9, zorder=2)
        ax.plot([radii[1],radii[2]],[zs[1],zs[2]], color=TH["link2_b"], lw=4, solid_capstyle='round', alpha=0.4, zorder=3)

       
        if abs(zs[2]-zs[3])>1:
            ax.plot([radii[2],radii[3]],[zs[2],zs[3]], color=TH["j3"], lw=4, ls='--', alpha=0.8, zorder=2)
            ax.annotate(f'd₃={self.d3.get():.0f}', (radii[2]+10,(zs[2]+zs[3])/2),
                        fontsize=7, color=TH["j3"],
                        bbox=dict(boxstyle='round,pad=0.2', facecolor=TH["card"],
                                  edgecolor=TH["j3"], alpha=0.8))

        
        jc=[TH["base_c"],TH["j1"],TH["j2"],TH["j3"],TH["tcp_c"]]
        jl=["Base","J1","J2","J3","TCP"]; jsz=[160,140,120,90,150]
        for i,(r,z) in enumerate(zip(radii,zs)):
            ax.scatter(r,z, color=jc[i], s=jsz[i]*2, alpha=0.08, zorder=4, edgecolors='none')
            ax.scatter(r,z, color=jc[i], s=jsz[i], zorder=5, edgecolors='white', linewidths=1.5)
            ax.annotate(jl[i], (r,z), textcoords="offset points", xytext=(10,10),
                        fontsize=7, color=jc[i], fontweight='bold', alpha=0.9)

        
        tcp_r = radii[4]; tcp_z = zs[4]
        GL = 10   
        GW = 20
        
        ax.plot([tcp_r-GW, tcp_r+GW], [tcp_z, tcp_z],
                color=TH["grip_c"], lw=5, solid_capstyle='round', alpha=0.95, zorder=4)
        
        for sign in [-1, +1]:
            rx = tcp_r + sign*GW
            ax.plot([rx, rx], [tcp_z, tcp_z-GL],
                    color=TH["grip_c"], lw=4, solid_capstyle='round', alpha=0.9, zorder=4)
           
            ax.scatter(rx, tcp_z-GL, color=TH["amber"], s=50, zorder=6,
                       edgecolors='white', linewidths=1.0)
        
        ax.annotate('', xy=(tcp_r+GW+4, tcp_z-GL),
                    xytext=(tcp_r-GW-4, tcp_z-GL),
                    arrowprops=dict(arrowstyle='<->', color=TH["grip_c"],
                                   lw=1.2, alpha=0.6))
        ax.text(tcp_r, tcp_z-GL-10, f"↕ {GL}mm", fontsize=6,
                color=TH["grip_c"], ha='center', alpha=0.7)

        if self.trace_on.get() and len(self.trace_pts)>1:
            tp=np.array(self.trace_pts)
            tr=np.sqrt(tp[:,0]**2+tp[:,1]**2); n=len(tp)
            for i in range(1,n):
                ax.plot(tr[i-1:i+1],tp[i-1:i+1,2],
                        color=TH["tcp_c"], alpha=0.05+0.8*(i/n), lw=1.5)

        ax.set_title("SCARA — Side View (RZ)", color=TH["text2"], fontsize=9, pad=6)

    
    def _reset(self):
        self.t1.set(0); self.t2.set(0); self.d3.set(0); self.t4.set(0)
        self._clear_trace()

    def _clear_trace(self):
        self.trace_pts.clear(); self._sched()

    def _toggle_anim(self):
        if self.animating:
            self.animating=False; self.btn_an.config(text="▶  Demo", fg=TH["cyan"])
        else:
            self.animating=True; self.btn_an.config(text="⏹  Stop", fg=TH["rose"])
            self.trace_on.set(True); self._anim(0)

    def _anim(self, s):
        if not self.animating: return
        t=s*0.03
        self.t1.set(round(60*np.sin(t),1))
        self.t2.set(round(45*np.sin(t*1.5+1),1))
        self.d3.set(round(100+80*np.sin(t*0.8),1))
        self.t4.set(round(90*np.sin(t*2),1))
        self.root.after(50, self._anim, s+1)

    def _fill_ik(self):
        try:
            fk=forward_kinematics(self.t1.get(),self.t2.get(),self.d3.get(),
                                   self.t4.get(),self.l1.get(),self.l2.get())
        except: return
        self.ik_x.set(round(fk['x'],2)); self.ik_y.set(round(fk['y'],2))
        self.ik_z.set(round(fk['z'],2)); self.ik_phi.set(round(fk['phi'],2))
        self.ik_st.config(text="✅ Đã lấy vị trí hiện tại.", fg=TH["green"])

    def _calc_ik(self):
        try:
            r=inverse_kinematics(self.ik_x.get(),self.ik_y.get(),self.ik_z.get(),
                                  self.ik_phi.get(),self.l1.get(),self.l2.get())
        except:
            self.ik_st.config(text="❌ Lỗi!", fg=TH["rose"]); return
        for w in self.ik_sf.winfo_children(): w.destroy()
        if "error" in r:
            self.ik_st.config(text=f"❌ {r['error']}", fg=TH["rose"]); self.ik_sols={}; return
        self.ik_sols=r
        self.ik_st.config(text="✅ Tìm thấy 2 nghiệm!", fg=TH["green"])
        for cn,lt,col in [("elbow_up","⬆ Elbow Up",TH["green"]),
                           ("elbow_down","⬇ Elbow Down",TH["purple"])]:
            sol=r[cn]
            cd=tk.Frame(self.ik_sf, bg=TH["input"], highlightthickness=1,
                        highlightbackground=TH["border"], padx=8, pady=4)
            cd.pack(fill="x", pady=2)
            tk.Label(cd, text=lt, font=self.F["badge"], bg=TH["input"], fg=col).pack(anchor="w")
            vf=tk.Frame(cd, bg=TH["input"]); vf.pack(fill="x", pady=(2,0))
            for param,val,u in [("θ₁",sol['t1'],"°"),("θ₂",sol['t2'],"°"),
                                 ("d₃",sol['d3'],"mm"),("θ₄",sol['t4'],"°")]:
                cell=tk.Frame(vf, bg=TH["input"]); cell.pack(side="left", expand=True)
                tk.Label(cell, text=param, font=self.F["tiny"], bg=TH["input"],
                         fg=TH["text3"]).pack()
                tk.Label(cell, text=f"{val:.1f}{u}", font=self.F["mono"],
                         bg=TH["input"], fg=TH["text"]).pack()

    def _apply_ik(self):
        cfg=self.ik_cfg.get()
        if not self.ik_sols:
            self.ik_st.config(text="⚠ Nhấn 'Tính IK' trước!", fg=TH["amber"]); return
        if cfg not in self.ik_sols:
            self.ik_st.config(text="❌ Không có nghiệm!", fg=TH["rose"]); return
        s=self.ik_sols[cfg]
        self.t1.set(s['t1']); self.t2.set(s['t2'])
        self.d3.set(s['d3']); self.t4.set(s['t4'])
        n="Elbow Up" if cfg=="elbow_up" else "Elbow Down"
        self.ik_st.config(text=f"✅ Đã áp dụng {n}!", fg=TH["green"])

   
    # TRAJECTORY
    
    def _traj_type_changed(self):
        if self.traj_type.get() == "line":
            self.traj_circ_f.pack_forget()
            self.traj_line_f.pack(fill="x")
        else:
            self.traj_line_f.pack_forget()
            self.traj_circ_f.pack(fill="x")

    def _run_trajectory(self):
        if self.traj_running:
            self.traj_running = False
            self.btn_traj.config(text="▶  Chạy quỹ đạo", bg="#059669")
            return

        l1v, l2v = self.l1.get(), self.l2.get()
        N = self.traj_N.get()
        T = self.traj_T.get()
        phi = self.traj_phi.get()

        
        if self.traj_type.get() == "line":
            pts = gen_line_trajectory(
                [self.traj_lx1.get(), self.traj_ly1.get(), self.traj_lz.get()],
                [self.traj_lx2.get(), self.traj_ly2.get(), self.traj_lz.get()],
                N=N, T=T)
        else:
            pts = gen_circle_trajectory(
                [self.traj_cx.get(), self.traj_cy.get()],
                self.traj_cr.get(), z=self.traj_cz.get(),
                N=N, T=T)

        
        joints = []
        for i, p in enumerate(pts):
            r = inverse_kinematics(p[0], p[1], p[2], phi, l1v, l2v)
            if "error" in r:
                self.traj_st.config(
                    text=f"❌ IK lỗi tại điểm {i}: {r['error']}", fg=TH["rose"])
                return
            sol = r["elbow_up"]
            joints.append([sol['t1'], sol['t2'], sol['d3'], sol['t4']])
        joints = np.array(joints)

        self.traj_pts_planned = pts
        self.traj_joints_log = joints
        self.traj_running = True
        self.trace_on.set(True)
        self._clear_trace()
        self.btn_traj.config(text="⏹  Dừng", bg=TH["rose"])
        self.traj_st.config(text=f"🏃 Đang chạy... {N} điểm, T={T}s", fg=TH["cyan"])

        dt = max(10, int(T * 1000 / N))  
        self._traj_step(joints, 0, dt)

    def _traj_step(self, joints, idx, dt):
        if not self.traj_running or idx >= len(joints):
            self.traj_running = False
            self.btn_traj.config(text="▶  Chạy quỹ đạo", bg="#059669")
            if idx >= len(joints):
                self.traj_st.config(text="✅ Hoàn thành quỹ đạo!", fg=TH["green"])
            return
        j = joints[idx]
        self.t1.set(round(j[0], 2))
        self.t2.set(round(j[1], 2))
        self.d3.set(round(j[2], 2))
        self.t4.set(round(j[3], 2))
        self.root.after(dt, self._traj_step, joints, idx+1, dt)

    def _show_traj_plots(self):
        """Show trajectory analysis plots in a new window."""
        if self.traj_joints_log is None or len(self.traj_joints_log) == 0:
            self.traj_st.config(text="⚠ Chạy quỹ đạo trước!", fg=TH["amber"])
            return
        joints = np.array(self.traj_joints_log)
        pts = self.traj_pts_planned
        T = self.traj_T.get()
        N = len(joints)
        t = np.linspace(0, T, N)

        win = tk.Toplevel(self.root)
        win.title("Phân Tích Quỹ Đạo — SCARA")
        win.configure(bg=TH["bg"])
        win.geometry("1000x700")

        fig = Figure(facecolor=TH["bg"], dpi=100)

        
        ax1 = fig.add_subplot(231, facecolor=TH["card"])
        ax1.plot(pts[:,0], pts[:,1], color=TH["cyan"], lw=2)
        ax1.plot(pts[0,0], pts[0,1], 'o', color=TH["green"], ms=8, label="Start")
        ax1.plot(pts[-1,0], pts[-1,1], 's', color=TH["rose"], ms=8, label="End")
        ax1.set_xlabel("X (mm)", color=TH["text3"], fontsize=8)
        ax1.set_ylabel("Y (mm)", color=TH["text3"], fontsize=8)
        ax1.set_title("Quỹ đạo XY", color=TH["text"], fontsize=10)
        ax1.set_aspect('equal')
        ax1.legend(fontsize=7)
        ax1.tick_params(colors=TH["text3"], labelsize=7)
        ax1.grid(True, color=TH["grid_c"], alpha=0.5)

        
        ax2 = fig.add_subplot(232, facecolor=TH["card"])
        ax2.plot(t, joints[:,0], color=TH["j1"], lw=1.5, label="θ₁")
        ax2.plot(t, joints[:,1], color=TH["j2"], lw=1.5, label="θ₂")
        ax2.plot(t, joints[:,3], color=TH["j4"], lw=1.5, label="θ₄")
        ax2.set_xlabel("t (s)", color=TH["text3"], fontsize=8)
        ax2.set_ylabel("Góc (°)", color=TH["text3"], fontsize=8)
        ax2.set_title("Góc khớp theo thời gian", color=TH["text"], fontsize=10)
        ax2.legend(fontsize=7)
        ax2.tick_params(colors=TH["text3"], labelsize=7)
        ax2.grid(True, color=TH["grid_c"], alpha=0.5)

        
        ax3 = fig.add_subplot(233, facecolor=TH["card"])
        ax3.plot(t, joints[:,2], color=TH["j3"], lw=1.5, label="d₃")
        ax3.set_xlabel("t (s)", color=TH["text3"], fontsize=8)
        ax3.set_ylabel("d₃ (mm)", color=TH["text3"], fontsize=8)
        ax3.set_title("d₃ theo thời gian", color=TH["text"], fontsize=10)
        ax3.legend(fontsize=7)
        ax3.tick_params(colors=TH["text3"], labelsize=7)
        ax3.grid(True, color=TH["grid_c"], alpha=0.5)

        
        dt_val = T / (N - 1)
        dq = np.diff(joints, axis=0) / dt_val
        tv = t[:-1] + dt_val/2
        ax4 = fig.add_subplot(234, facecolor=TH["card"])
        ax4.plot(tv, dq[:,0], color=TH["j1"], lw=1.5, label="ω₁")
        ax4.plot(tv, dq[:,1], color=TH["j2"], lw=1.5, label="ω₂")
        ax4.plot(tv, dq[:,3], color=TH["j4"], lw=1.5, label="ω₄")
        ax4.set_xlabel("t (s)", color=TH["text3"], fontsize=8)
        ax4.set_ylabel("Vận tốc (°/s)", color=TH["text3"], fontsize=8)
        ax4.set_title("Vận tốc khớp", color=TH["text"], fontsize=10)
        ax4.legend(fontsize=7)
        ax4.tick_params(colors=TH["text3"], labelsize=7)
        ax4.grid(True, color=TH["grid_c"], alpha=0.5)

        
        dp = np.diff(pts, axis=0) / dt_val
        v_tcp = np.sqrt(dp[:,0]**2 + dp[:,1]**2 + dp[:,2]**2)
        ax5 = fig.add_subplot(235, facecolor=TH["card"])
        ax5.plot(tv, v_tcp, color=TH["tcp_c"], lw=2)
        ax5.set_xlabel("t (s)", color=TH["text3"], fontsize=8)
        ax5.set_ylabel("v (mm/s)", color=TH["text3"], fontsize=8)
        ax5.set_title("Vận tốc TCP (task space)", color=TH["text"], fontsize=10)
        ax5.tick_params(colors=TH["text3"], labelsize=7)
        ax5.grid(True, color=TH["grid_c"], alpha=0.5)

        
        l1v, l2v = self.l1.get(), self.l2.get()
        err = []
        for i in range(N):
            fk = forward_kinematics(joints[i,0], joints[i,1], joints[i,2], joints[i,3], l1v, l2v)
            ex = fk['x'] - pts[i,0]
            ey = fk['y'] - pts[i,1]
            ez = fk['z'] - pts[i,2]
            err.append(np.sqrt(ex**2 + ey**2 + ez**2))
        ax6 = fig.add_subplot(236, facecolor=TH["card"])
        ax6.plot(t, err, color=TH["amber"], lw=1.5)
        ax6.set_xlabel("t (s)", color=TH["text3"], fontsize=8)
        ax6.set_ylabel("Sai số (mm)", color=TH["text3"], fontsize=8)
        ax6.set_title("Sai số IK tracking", color=TH["text"], fontsize=10)
        ax6.tick_params(colors=TH["text3"], labelsize=7)
        ax6.grid(True, color=TH["grid_c"], alpha=0.5)

        fig.tight_layout(pad=2)
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)



if __name__ == "__main__":
    root = tk.Tk()
    try:
        import ctypes; root.update()
        h=ctypes.windll.user32.GetParent(root.winfo_id())
        v=ctypes.c_int(2)
        ctypes.windll.dwmapi.DwmSetWindowAttribute(h,20,ctypes.byref(v),ctypes.sizeof(v))
    except: pass
    App(root)
    root.mainloop()
