
# try this code , i think you have to do some modifications too 
R3_6 = R0_3.transpose() * R0_EE
x4 = -R3_6[0,2]
y6 = -R3_6[1,1]
y5 = sqrt(R3_6[0,2]*R3_6[0,2] + R3_6[2,2]*R3_6[2,2])

theta4 = atan2(R3_6[2,2], x4)
theta4 = theta4.evalf()

# theta5 = acos(R3_6[1,2])
theta5 = atan2(y5,R3_6[1,2])
theta5 = theta5.evalf()

theta6 = atan2(y6, R3_6[1,0])
theta6 = theta6.evalf()
