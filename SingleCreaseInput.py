"""
Created on June 21 2022
@author: Zhang Qian
"""


from math import *
import numpy as np

class SingleCrease(object):
    """
    A single crease class
    This procedure is used for the single crease simulation for Edge Displacement actuator.
    Tension and Compression.
    Crease Stiffness
            *          ---
           * *         ---
          *   *        ---
         *     *       ---
        *       *      ---
    D<-*         *->D  ---
    There are two plates connected by one creases.
    """
    def __init__(self, jobname, Stiffness, dx, width, height, SectorAngle, InitialAngle, thickness, foldingstate, Loadregion):
        self.jobname = jobname
        self.stiffness = Stiffness
        self.dx = dx
        self.width = width
        self.height = height
        self.SectorAngle = SectorAngle
        self.InitialAngle = InitialAngle
        self.thickness = thickness
        self.foldingstate = foldingstate
        self.loadregion = Loadregion


    def _coordinate(self, geometry, InitialAngle, PlateNumber):
        """
        The vertex number of the plate starts from the left-up node,
        then increases count-clockwise, like
        1, 4
        2, 3
        The node 3 of left plate is selected as the origin point.
        the line 3/4 of left plate is the y-axis
        the 3/4 of left is connected by 2/1 of right plate
        :arg geometry [width, Height and sector angle] of the plate.
        :arg float InitialAngle: the initial angle of the two plate.
        :arg int PlateNumber: 0 for left plate and 1 for right plate
        """
        Rotation = (pi - InitialAngle)/2
        T = np.array([[cos(Rotation), 0, sin(Rotation)],
                      [0, 1, 0],
                      [-sin(Rotation), 0, cos(Rotation)]])
        width = geometry[0]
        height = geometry[1]
        SecAngle = geometry[2]
        CoordR = np.ones((4, 3))
        CoordR[0] = np.array([0, height, 0])
        CoordR[1] = np.array([0, 0, 0])
        temp3 = np.array([width * sin(SecAngle), -width * cos(SecAngle), 0])
        CoordR[2] = np.dot(T, temp3)
        CoordR[3] = CoordR[2] + CoordR[0]
        if PlateNumber == 1:
            Coord = CoordR
        else:
            CoordL = np.ones((4, 3))
            M = np.array([[-1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
            CoordL[0] = np.dot(M, CoordR[3])
            CoordL[1] = np.dot(M, CoordR[2])
            CoordL[2] = CoordR[1]
            CoordL[3] = CoordR[0]
            Coord = CoordL

        return Coord


    def _Mesh(self, nw, nh, Coord, identifier):
        """
        Only support Quadrilateral mesh
        Calculate the node list and element list in global system
        nw: width/dx
        nh: height/dx
        identifier: the plate number, int
        :return: NodeCoord Node list A(,4), (Node Number, x, y, z)
        :return: Element element lsit B(,5),
                (Element Number, Node1, Node2, Node 3, Node4 )
        :return: EdgeSet --Edge Node set for Connection design.
        Edge Number: 0 is for line 1-2, 1 is for line 2-3. 2 is for line 4-3
                3 is for line 1-4. The order of the line is considered.
        """
        NodeCoord = np.zeros(((nw + 1) * (nh + 1), 4))
        Element = np.zeros((nw * nh, 5), dtype=np.int32)
        EdgeSet = []
        Number = 1
        for i in range(nw + 1):
            start = Coord[0] + (i / nw) * (Coord[3] - Coord[0])
            end = Coord[1] + (i / nw) * (Coord[2] - Coord[1])
            for j in range(nh + 1):
                temp = start + (j / nh) * (end - start)
                Global_Number = Number + identifier * (nw + 1) * (nh + 1)
                NodeCoord[Number - 1] = np.array([Global_Number, temp[0], temp[1], temp[2]])
                Number += 1

        Number = 1
        for i in range(nw):
            for j in range(nh):
                Global_symbol = identifier * (nw + 1) * (nh + 1)
                Element[Number - 1] = np.array([Number + identifier * nw * nh,
                                                Global_symbol + i * (nh + 1) + j + 1,
                                                Global_symbol + i * (nh + 1) + j + 2,
                                                Global_symbol + (i + 1) * (nh + 1) + j + 2,
                                                Global_symbol + (i + 1) * (nh + 1) + j + 1])
                Number += 1
        for i in range(4):
            plateIncrease = identifier * (nw + 1) * (nh + 1)
            if (i % 2) == 0:
                EdgeIncrease = 1 if i == 0 else nw * (nh + 1) + 1
                EdgeNode = [itme for itme in range(plateIncrease + EdgeIncrease, plateIncrease + EdgeIncrease + nh + 1)]
                EdgeSet.append(EdgeNode)
            else:
                EdgeIncrease = nh + 1 if i == 1 else 1
                EdgeNode = [itme for itme in range(plateIncrease + EdgeIncrease, plateIncrease + EdgeIncrease + (nh + 1) * nw + 1, nh + 1)]
                EdgeSet.append(EdgeNode)

        return NodeCoord, Element, EdgeSet

    def Inpwrite(self):
        """
        Basic Information Input Module
        The plate geometry, [width, height, Sector Angle]
        the crease length is the height, and  Sector angle between line 1/4 to 1/2, default = 90
        InitialAngle: Define the initial rest angle of the single crease
        For single crease, there are 2 plates.
        """
        width = self.width
        height = self.height
        SectorAngle = self.SectorAngle
        geometry = np.array([width, height, SectorAngle])
        InitialAngle = self.InitialAngle
        PlateNumber = 2

        """
        Mesh Information
        dx is the size of mesh, nw, nh are the number of seed of
        edge with the width and height.
        NodeNumber and ElementNumber are the number of the nodes and 
        elements for each plate.
        dx > 10 * thickness for shell element
        """
        dx = self.dx
        nw = floor(width / dx)
        nh = floor(height / dx)
        NodeNumber = (nw + 1) * (nh + 1)
        ElementNumber = nw * nh


        """
        Node List and Element List Modulus
        Assembly.
        Including NodeCoord and Element Information for all plate.
        The geometry points for each plate is defined as
        1  4
        2  3
        Plate Number: 0 is for left plate, while 1 is for right plate
        Edge Number: 0 is for line 1-2, 1 is for line 2-3. 2 is for line 4-3 
                    3 is for line 1-4. The order of the line is considered.
        Position of the edge node in the mesh, Plate Number, Edge Number, Node Number 
        such as Plate Number = 0, Edge Number = 1, Node Number = 5 means the nodes is located
        at the left plate, line 2-3, the fifth node. 
        """
        NodeCoord = []
        Element = []
        EdgeSet = []
        for i in range(PlateNumber):
            Coord_temp = self._coordinate(geometry, InitialAngle, i)
            NodeCoord_temp, Element_temp, Edge_temp = self._Mesh(nw, nh, Coord_temp, i)
            if i == 0:
                NodeCoord = NodeCoord_temp
                Element = Element_temp
                EdgeSet = Edge_temp
            else:
                NodeCoord = np.concatenate((NodeCoord, NodeCoord_temp))
                Element = np.concatenate((Element, Element_temp))
                EdgeSet = EdgeSet + Edge_temp


        """
        Material Information of FEM Model
        The plate is elastic material
        Including Thickness, Elastic Module, Poisson ratio 
        unit [mm], [kg] [N] [rad]
        PET material except for the thickness
        ref: Embedded Actuation for Shape-Adaptive Origami, 2021, k=0.05
        B = Et^3/12(1-v^2)
        """
        thickness = self.thickness
        Elastic_Module = 3200.0
        Possion_ratio = 0.43
        Stiffness = self.stiffness
        CreaseStiffness = Stiffness * dx
        B = Elastic_Module * thickness ** 3 / 12 / (1 - Possion_ratio ** 2)
        ratio = B / Stiffness
        print("The length scale B/k is ", ratio)
        Total_NodeNumber = PlateNumber * NodeNumber
        Total_ElementNumber = PlateNumber * ElementNumber

        """
        Loading Information
        According to the edge number in the Total system,
        Transform: global = 4 * plate + local 
        For example (1,2) corresponds to the the second plate and the third edge.
        The number of the global system is 4 * 1 + 2 = 6
        In self-folding:
        Load_Edge_left U1 U2 U3 are constrained
        Load_Edge_right U2 and U3 are constrained 
        """
        Load_Edge_left = 4 * 0 + 0
        Load_Edge_right = 4 * 1 + 2

        # print(EdgeSet[Load_Edge_left])
        lowerbound = ceil((0.5-self.loadregion/2.0) * len(EdgeSet[Load_Edge_left]))
        uplowerbound = floor((0.5+self.loadregion/2.0) * len(EdgeSet[Load_Edge_left]))
        LoadLeft = EdgeSet[Load_Edge_left][lowerbound:uplowerbound]
        LoadRight = EdgeSet[Load_Edge_right][lowerbound:uplowerbound]
        # print(LoadLeft)
        # # print(LoadRight)


        if self.foldingstate == 1:
            Dis = width * sin(SectorAngle) * (1 - sin(InitialAngle / 2))
        elif self.foldingstate == 0:
            Dis = -width * sin(SectorAngle) * (sin(InitialAngle / 2))
        else:
            print("error for the folding state input")

        """
        Connection Information
        Use the specified method instead of the query method
        For the two plate, the connection Information can be given 
        as (0,2) to (1,0) between the third edge of first plate and 
        the first edge of the second plate
        The local cooridnate system is also defined through the node information
        There is only one connection edge
        """
        ConnectionNumber = 1
        PlatePair = [0, 1]
        EdgePair = [2, 0]
        EdgePair_global = [4 * x + y for x, y in zip(PlatePair, EdgePair)]
        print("The Connection Pair (Edge Number) is ", EdgePair_global)


        """
        Output Information for analysis
        Configuration Label
        OutputX : Plate0 the middle line along the width direction.
        OutputY : Plate0 the connection edge (Number: 2)
        """
        Label = EdgePair_global[0]
        OutputY = EdgeSet[Label]
        # print(nh)
        OutputX = []
        for x in range(nw + 1):
            OutputX.append(ceil((nh + 1) / 2) + x * (nh + 1))
        with open('OutputX.txt', "w") as f:
            for i in range(len(OutputX)):
                f.write(str(OutputX[i]) + "\n")
        with open('OutputY.txt', "w") as f:
            for i in range(len(OutputY)):
                f.write(str(OutputY[i]) + "\n")
        """
        Inp file Module
        Label information: Width, Height, Sector Angle, Initial Angle, PlateNumber
        Define the name of the inp file as "singleCrease"
        Static calculation in General/static solver. 
        """

        inp_file = open(self.jobname+".inp", "w")
        inp_file.write("*Heading\n")
        inp_file.write("**Job Name and Model Name:SingleCrease\n")
        inp_file.write("*Preprint, echo=NO, model=NO, history=NO, contact=NO\n")
        inp_file.write("**\n")

        inp_file.write("**PARTS\n")
        inp_file.write("*Part,name=Crease\n")
        inp_file.write("*Node\n")
        # for the node list of left plate
        for i in range(Total_NodeNumber):
            inp_file.write(str(i + 1) + "," + str(NodeCoord[i][1]) + ","
                           + str(NodeCoord[i][2]) +"," + str(NodeCoord[i][3]) +"\n")
        inp_file.write("**\n")

        inp_file.write("*Element,type=S4R\n")
        for i in range(Total_ElementNumber):
            inp_file.write(str(i + 1) + "," + str(Element[i][1]) + "," + str(Element[i][2]) + ","
                           + str(Element[i][3]) +"," + str(Element[i][4]) +"\n")
        inp_file.write("**\n")

        inp_file.write("*Nset, nset=SET-All, generate\n")
        inp_file.write(" 1," + str(Total_NodeNumber) + ",1\n")
        inp_file.write("*Elset, elset=SET-All, generate\n")
        inp_file.write("1," + str(Total_ElementNumber) + ",1\n")

        inp_file.write("** Section: \n")
        inp_file.write("*Shell Section, elset=SET-All, material=Self-define\n")
        inp_file.write(str(thickness)+", 5\n")
        inp_file.write("*End Part\n")

        """
        Assembly Modulus, including the node set definition
        element set definition, and connection definition.
        """

        inp_file.write("** ASSEMBLY\n")
        inp_file.write("*Assembly, name=Assembly\n")
        inp_file.write("*Instance, name=SingleCreaseModel, part=Crease\n")
        inp_file.write("*End Instance\n")

        inp_file.write("*Nset, nset=Set-Left, instance=SingleCreaseModel\n")
        length = len(LoadLeft)
        for item in range(length - 1):
            if item % 10 == 9:
                inp_file.write("\n")
            inp_file.write(str(LoadLeft[item])+",")
        inp_file.write(str(LoadLeft[-1])+"\n")

        inp_file.write("*Nset, nset=Set-Right, instance=SingleCreaseModel\n")
        length = len(LoadRight)
        for item in range(length - 1):
            if item % 10 == 9:
                inp_file.write("\n")
            inp_file.write(str(LoadRight[item])+",")
        inp_file.write(str(LoadRight[-1])+"\n")

        inp_file.write("*Elset, elset=SET-Plate, instance=SingleCreaseModel, generate\n")
        inp_file.write("1," + str(Total_ElementNumber) + ",1\n")

        # there is only one connection Pair in this analysis
        inp_file.write("*Element, type=CONN3D2\n")
        Pair1 = EdgePair_global[0]
        Pair2 = EdgePair_global[1]
        PairLen = len(EdgeSet[Pair1])
        for i in range(PairLen):
            inp_file.write(str(i+1) + ", SingleCreaseModel." + str(EdgeSet[Pair1][i]) + ", SingleCreaseModel." + str(EdgeSet[Pair2][i]) + "\n")
        inp_file.write("*Nset, nset=Set-Connection, instance=SingleCreaseModel\n")
        for item in range(PairLen-1):
            if item % 6 == 5:
                inp_file.write("\n")
            inp_file.write(str(EdgeSet[Pair1][item])+","+str(EdgeSet[Pair2][item])+",")
        inp_file.write(str(EdgeSet[Pair1][-1])+","+str(EdgeSet[Pair2][-1])+"\n")
        inp_file.write("*Elset, elset=Set-Connection, generate\n")
        inp_file.write("1," + str(PairLen) + ",1 \n")
        # inp_file.write("*Orientation, name=csys-Con, DEFINITION=NODES\n")
        # inp_file.write("SingleCreaseModel." + str(EdgeSet[Pair1][0]) + ", SingleCreaseModel.1, SingleCreaseModel." + str(EdgeSet[Pair1][-1]) + "\n")
        inp_file.write("*Orientation, name=csys-Con\n")
        inp_file.write("0.,1.0,0.,-1.0,0.,0.\n")
        inp_file.write("1,0.\n")
        inp_file.write("*Connector Section, elset=Set-Connection, behavior=ConnProp-1\n")
        inp_file.write("Hinge,\n")
        inp_file.write("csys-Con,\n")
        inp_file.write("*End Assembly\n")
        inp_file.write("*Connector Behavior, name=ConnProp-1\n")
        inp_file.write("*Connector Elasticity, component=4\n")
        inp_file.write(str(CreaseStiffness)+",\n")

        """
        Standard Module
        Including Material Information, Step Information, 
        Boundary Information
        """
        inp_file.write("** MATERIALS\n")
        inp_file.write("*Material, name=Self-define\n")
        inp_file.write("*Elastic\n")
        inp_file.write(str(Elastic_Module) + "," + str(Possion_ratio) + "\n")

        inp_file.write("** STEP: Step-1\n")
        inp_file.write("*Step, name=Step-1, nlgeom=YES, inc=1000\n")
        inp_file.write("*Static\n")
        inp_file.write("0.01, 1., 1e-15, 0.01\n")

        inp_file.write("** BOUNDARY CONDITIONS\n")
        inp_file.write("** Name: BC-LEFT Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-left, 1, 1, " + str(-Dis) + "\n")
        inp_file.write("Set-left, 2, 2\n")
        inp_file.write("Set-left, 3, 3\n")
        inp_file.write("** Name: BC-RIGHT Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-right, 1, 1, " + str(Dis) + "\n")
        inp_file.write("Set-right, 2, 2\n")
        inp_file.write("Set-right, 3, 3\n")

        inp_file.write("** CONTROLS\n")
        inp_file.write("*Controls, reset\n")
        inp_file.write("*Controls, parameters=time incrementation\n")
        inp_file.write("8, 10, , , , , , 50, , , \n")
        inp_file.write("** OUTPUT REQUESTS\n")
        inp_file.write("*Restart, write, frequency=0\n")
        inp_file.write("** FIELD OUTPUT: F-Output-1\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Node Output\n")
        inp_file.write("CF, COORD, RF, U\n")
        inp_file.write("*Element Output, directions=YES\n")
        inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        inp_file.write("** FIELD OUTPUT: F-Output-2\n")
        inp_file.write("**\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Element Output, elset=Set-Connection, directions=YES\n")
        inp_file.write("CTF,\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-1\n")
        inp_file.write("*Output, history, variable=PRESELECT\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-2\n")
        inp_file.write("*Output, history\n")
        inp_file.write("*Element Output, elset=Set-Connection\n")
        inp_file.write("CTM1, \n")
        inp_file.write("** HISTORY OUTPUT: H-Output-3\n")
        inp_file.write("*Output, history\n")
        inp_file.write("*Energy Output, elset=SET-Plate\n")
        inp_file.write("ALLSE, \n")
        inp_file.write("*End Step\n")
        #
        # inp_file.write("** STEP: Step-2\n")
        # inp_file.write("*Step, name=Step-2, nlgeom=YES, inc=1000\n")
        # inp_file.write("*Static\n")
        # inp_file.write("0.01, 1., 1e-15, 0.01\n")
        #
        # inp_file.write("** BOUNDARY CONDITIONS\n")
        # inp_file.write("** Name: BC-LEFT Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-left, 1, 1, " + str(-Dis * 1.1) + "\n")
        # inp_file.write("** Name: BC-RIGHT Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-right, 1, 1, " + str(Dis * 1.1) + "\n")
        #
        # inp_file.write("** CONTROLS\n")
        # inp_file.write("*Controls, reset\n")
        # inp_file.write("*Controls, parameters=time incrementation\n")
        # inp_file.write("8, 10, , , , , , 50, , , \n")
        # inp_file.write("** OUTPUT REQUESTS\n")
        # inp_file.write("*Restart, write, frequency=0\n")
        # inp_file.write("** FIELD OUTPUT: F-Output-1\n")
        # inp_file.write("*Output, field\n")
        # inp_file.write("*Node Output\n")
        # inp_file.write("CF, COORD, RF, U\n")
        # inp_file.write("*Element Output, directions=YES\n")
        # inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        # inp_file.write("** FIELD OUTPUT: F-Output-2\n")
        # inp_file.write("**\n")
        # inp_file.write("*Output, field\n")
        # inp_file.write("*Element Output, elset=Set-Connection, directions=YES\n")
        # inp_file.write("CTF,\n")
        # inp_file.write("** HISTORY OUTPUT: H-Output-1\n")
        # inp_file.write("*Output, history, variable=PRESELECT\n")
        # inp_file.write("** HISTORY OUTPUT: H-Output-2\n")
        # inp_file.write("*Output, history\n")
        # inp_file.write("*Element Output, elset=Set-Connection\n")
        # inp_file.write("CTM1, \n")
        # inp_file.write("** HISTORY OUTPUT: H-Output-3\n")
        # inp_file.write("*Output, history\n")
        # inp_file.write("*Energy Output, elset=SET-Plate\n")
        # inp_file.write("ALLSE, \n")
        # inp_file.write("*End Step\n")
        inp_file.close()
