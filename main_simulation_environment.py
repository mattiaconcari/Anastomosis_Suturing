#%% buca o non buca
import pybullet as p
import pybullet_data
import numpy as np
import time
from pathlib import Path

# Connessione a PyBullet
p.connect(p.GUI)
p.setAdditionalSearchPath(pybullet_data.getDataPath())
p.resetSimulation(p.RESET_USE_DEFORMABLE_WORLD)
p.setGravity(0, 0, -9.81)

BASE_DIR = Path(__file__).resolve().parent  # Cartella dello script
ASSETS_DIR = BASE_DIR / "assets"  # Cartella contenente i file

# File di input
STL_FILE = str(ASSETS_DIR / "ago_66mm.stl")  # Converti in stringa per PyBullet
VESSEL_OBJ = str(ASSETS_DIR / "Arteria_piena.obj")
VESSEL_VTK = str(ASSETS_DIR / "Arteria_piena.vtk")
# Caricamento del piano
p.loadURDF("plane.urdf", [0, 0, 0])

# Funzione per caricare un oggetto da un file STL
def create_needle_from_stl(file_path, position, orientation, scale=0.01):
    collision_shape_id = p.createCollisionShape(shapeType=p.GEOM_MESH, fileName=file_path, meshScale=[scale, scale, scale])
    visual_shape_id = p.createVisualShape(shapeType=p.GEOM_MESH, fileName=file_path, meshScale=[scale, scale, scale])
    needle_id = p.createMultiBody(baseMass=0.1, baseCollisionShapeIndex=collision_shape_id,
                                  baseVisualShapeIndex=visual_shape_id, basePosition=position,
                                  baseOrientation=orientation)
    return needle_id


# Parametri del vaso
vessel_outer_radius = 0.5  # Raggio esterno del vaso (50 mm)
vessel_length = 1.0  # Lunghezza del vaso


# Parametri del semicerchio
radius = 0.33  # Raggio del semicerchio in metri (33 mm)
angle = 0  # Angolo iniziale
needle_center=[0, 0, 2*vessel_outer_radius]  # Posizione del centro fisso
# Parametri dell'ago
needle_position = [0, 0, 2*vessel_outer_radius]  # Posizione iniziale
needle_orientation = p.getQuaternionFromEuler([np.pi, 0, 0])  # Orientamento iniziale
needle_id = create_needle_from_stl(STL_FILE, needle_position, needle_orientation, scale=0.01)


# Posizionamento del vaso
vessel1 = p.loadSoftBody(
    VESSEL_OBJ,
    simFileName=VESSEL_VTK,
    basePosition=[-vessel_outer_radius, -0.25, 0],  # Centro del vaso allineato sull'asse Z
    baseOrientation=p.getQuaternionFromEuler([np.pi/2, 0, 0]),
    scale=0.1,
    mass=4,
    useNeoHookean=1,
    NeoHookeanMu=5000,
    NeoHookeanLambda=1000,
    NeoHookeanDamping=0.01,
    useSelfCollision=1,
    frictionCoeff=0.5,
    collisionMargin=0.01
)

"""vessel2=p.loadSoftBody(VESSEL_OBJ, simFileName=VESSEL_VTK, basePosition=[-vessel_outer_radius, 2*vessel_length+0.025, 0], baseOrientation=p.getQuaternionFromEuler([ np.pi/2,0, 0]),
                       scale=0.1, mass=4, useNeoHookean=1, NeoHookeanMu=5000, NeoHookeanLambda=1000, NeoHookeanDamping=0.01, useSelfCollision=1, frictionCoeff=0.5, collisionMargin=0.01)
"""

# Slider per controllare velocità di rotazione e configurazione della telecamera
rotation_speed_slider = p.addUserDebugParameter("Rotation Speed", -5.0, 10.0, 0.5)
camera_yaw_slider = p.addUserDebugParameter("Camera Yaw", 0, 360, 145)
camera_pitch_slider = p.addUserDebugParameter("Camera Pitch", -180, 180, -5)
camera_distance_slider = p.addUserDebugParameter("Camera Distance", 0.1, 10, 2)

# Centro del semicerchio allineato con il centro del vaso
semicircle_center = [0, 0, 2*vessel_outer_radius]  # Centro fisso del semicerchio

num_suture=6
current_suture=1
suture_angle=2*np.pi/num_suture
rot_angle=0
movement=0

# Ciclo di simulazione
while p.isConnected():
    # Simulazione
    p.stepSimulation()

    # Ottieni valori dagli sliders
    rotation_speed = p.readUserDebugParameter(rotation_speed_slider)
    yaw = p.readUserDebugParameter(camera_yaw_slider)
    pitch = p.readUserDebugParameter(camera_pitch_slider)
    distance = p.readUserDebugParameter(camera_distance_slider)


    # Aumenta l'angolo di rotazione
    angle += rotation_speed * (1 / 240)  # Incremento angolare
    if current_suture==1:
        needle_orientation=p.getQuaternionFromEuler([np.pi+angle,0,0])
        p.resetBasePositionAndOrientation(needle_id, needle_center, needle_orientation)
    
    while movement==1: #mi serve per far muovere l'ago e non farlo teletrasportare
        p.stepSimulation()
        angle += rotation_speed * (1 / 240)
        if angle <= current_angle:
            p.resetBasePositionAndOrientation(needle_id, last_needle_center, needle_orientation)
            needle_orientation=p.getQuaternionFromEuler([np.pi,current_angle,0])
            needle_center=[vessel_outer_radius*np.sin(previous_current_angle+angle),0,vessel_outer_radius+vessel_outer_radius*np.cos(previous_current_angle+angle)]#faccio luovere da dove ho lasciato
            p.resetBasePositionAndOrientation(needle_id, needle_center, needle_orientation)
        else:
            needle_center=[vessel_outer_radius*np.sin(current_angle),0,vessel_outer_radius+vessel_outer_radius*np.cos(current_angle)]
            movement=0 #finito il movimento si ferma anche questo pezzo della macchina a stati

    while rot_angle == 1:
        angle += rotation_speed * (1 / 240)
        needle_orientation=p.getQuaternionFromEuler([np.pi+angle,current_angle,0])
        needle_center=[vessel_outer_radius*np.sin(current_angle),0,vessel_outer_radius+vessel_outer_radius*np.cos(current_angle)]
        p.stepSimulation()
        p.resetBasePositionAndOrientation(needle_id, needle_center, needle_orientation)
        rotation_speed = p.readUserDebugParameter(rotation_speed_slider)
        yaw = p.readUserDebugParameter(camera_yaw_slider)
        pitch = p.readUserDebugParameter(camera_pitch_slider)
        distance = p.readUserDebugParameter(camera_distance_slider)
        p.resetDebugVisualizerCamera(
            cameraDistance=distance,
            cameraYaw=yaw,
            cameraPitch=pitch,
            cameraTargetPosition=semicircle_center
        )
        if angle>33/18*np.pi:
            rot_angle=0
        #  angle=np.pi+0.01

    if angle>33/18*np.pi and current_suture<=num_suture and rot_angle==0:
        angle=0
        previous_current_angle=suture_angle*(current_suture-1)#mi serve per non cominciare sempre da sopra il vaso ma da dove ave lasciato
        
        current_angle=current_suture*suture_angle
        #needle_orientation=p.getQuaternionFromEuler([np.pi+angle,current_angle,0])
        last_needle_center=[1+vessel_outer_radius*np.sin(previous_current_angle),0,1+vessel_outer_radius+vessel_outer_radius*np.cos(previous_current_angle)]
        current_suture=current_suture+1
        rot_angle=1
        movement=1

    if current_suture>num_suture:
        needle_position = [0, 0, 2*vessel_outer_radius]  # Posizione iniziale
        needle_orientation = p.getQuaternionFromEuler([0, 0, 0])  # Orientamento iniziale
        p.resetBasePositionAndOrientation(needle_id, needle_center, needle_orientation)
    
    #angle += rotation_speed * (1 / 240)
    p.resetBasePositionAndOrientation(needle_id, needle_center, needle_orientation)





    # Configura telecamera
    p.resetDebugVisualizerCamera(
        cameraDistance=distance,
        cameraYaw=yaw,
        cameraPitch=pitch,
        cameraTargetPosition=semicircle_center
    )

    time.sleep(1 / 240)

