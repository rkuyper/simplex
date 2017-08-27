#This code implements the simplex algorithm in object-oriented Python, using the Devex pivot rule
#Written to serve as a sample of my programming ability

import copy

#Class: Linear Inequality
#Attributes:
##a: contains the list defining the coefficients
##op: either "==", "<=" or ">=", defining which way the inequality goes
##b: contains the right-hand side of the inequality
##n: the number of variables
#Procedures:
##evaluate(x): evaluates the value of a . x
##
class LinIneq(object):

    def __init__(self, a, op, b):
        self.a = a
        self.op = op
        self.b = b
        self.n = len(a)

    #Procedure evaluate
    #input: x
    #output: the inner product of x and a
    #undefined values are assumed to be 0
    def evaluate(self,x):
        sum = 0
        for i in range(0,min(self.n,len(x))):
            sum += self.a[i]*x[i]
        if self.op == "==" and sum == self.b:
            return True
        if self.op == "<=" and sum <= self.b:
            return True
        if self.op == ">=" and sum >= self.b:
            return True
        else:
            return False

#Class: Linear System, a system of Linear Inequalities
#Attributes:
##eqns: a list of LinIneq
##m: the number of LinIneq
#Procedures:
##is_feasible(x): tells us whether x is a feasible solution,
##                i.e. whether it satisfies the inequalities in eqns
##equalise_length: updates the system of equations such that they are all
##                 of the same length, i.e. same number of variables
class LinSystem(object):

    def __init__(self,eqns):
        self.eqns = eqns
        self.m = len(eqns)

    #Procedure is_feasible(x)
    #input: a list x
    #output: a Bool telling us whether all the linear inequalities in self.eqns
    #are satisfied by x. Undefined values are assumed to be 0
    def is_feasible(self,x):
        for i in range(0,self.m):
            if not self.eqns[i].evaluate(x):
                return False
        return True

    #Procedure equalise_lengths
    #modifies eqns to all have the same length, namely the maximum of the existing lengths
    def equalise_lengths(self):
        new_n = max([eqn.n for eqn in self.eqns])

        for i in range(0,self.m):
            self.eqns[i].a += ([0]*(new_n-self.eqns[i].n))
            self.eqns[i].n = new_n

#Class Linear Program
#Attributes:
##system: an object of type LinSystem
##objective: a list giving the coefficients for the objective
##minmax: either "min" or "max", telling us whether we maximise of minimise
#Procedures:
##evaluate(x): evaluates the objective at x
##normalise(x): puts the linear program into normal form
##solve(x): solves the linear program
##printmatrix: prints the matrix representing the linear program, for debugging purposes

class LinProgram(object):

    def __init__(self,system,objective,minmax):
        self.system = system
        self.objective = objective
        self.minmax = minmax

    #Procedure evaluate(x)
    #input: a list
    #output: the objective evaluated at x, i.e. the inner product of self.objective and x
    #undefined values are assumed to be 0
    def evaluate(self,x):
        sum = 0
        for i in range(0,max(len(self.objective),len(x))):
            sum += self.objective[i]*x[i]
        return sum

    #Procedure normalise
    #Transforms the linear program into one of standard form, i.e.:
    ##-The target is to minimise
    ##-All inequalities are of the form >= 0
    ##-All variables have an inequality >= 0
    def normalise(self):
        #First, make sure everything has the same length
        self.system.equalise_lengths()

        n = self.system.eqns[0].n
        m = self.system.m

        for i in range(0,self.system.m):
            #Split every unrestricted variable into the difference of two positive variables
            self.system.eqns[i].a += [-x for x in self.system.eqns[i].a] + [0.0]*m
            self.system.eqns[i].n = 2*n+m

            thisa = [0.0]*(2*n+m)

            #Add a slack variable
            thisa[2*n+i] = 1
            self.system.eqns.append(LinIneq(thisa,">=",0))

            if self.system.eqns[i].op == "<=":
                self.system.eqns[i].a[2*n+i] = 1.0
                self.system.eqns[i].op = "=="
            if self.system.eqns[i].op == ">=":
                self.system.eqns[i].a = [-u for u in self.system.eqns[i].a]
                self.system.eqns[i].b *= -1.0
                self.system.eqns[i].a[2*n+i] = 1.0
                self.system.eqns[i].op = "=="

        self.system.m *= 2

        #Change maximisation into minimisation
        if self.minmax == "max":
            self.minmax = "min"
            self.objective = [-i for i in self.objective]

        #Update the objective
        self.objective = self.objective + [-x for x in self.objective] + [0]*m

    #Procedure solve
    #Solve the linear system using a two-phase simplex algorithm
    #Output: either the optimal solution, the string "Unbounded" or the string "No feasible solution"
    def solve(self):
        n = self.system.eqns[0].n
        m = self.system.m

        #Make a copy to preserve the original program
        cpy = copy.deepcopy(self)
        cpy.normalise()

        #Phase 1, finding a feasible solution. We add a variable for every inequality
        phase1 = copy.deepcopy(cpy)
        init_sol = [0.0]*phase1.system.eqns[0].n
        for i in range(0, phase1.system.m):
            phase1.system.eqns[i].a += [0]*phase1.system.m

        for i in range(0,phase1.system.m):
            #Make sure that we have +1 in basic columns
            if phase1.system.eqns[i].b < 0:
                phase1.system.eqns[i].a = [-x for x in phase1.system.eqns[i].a]
                phase1.system.eqns[i].b *= -1.0

            phase1.system.eqns[i].a[phase1.system.eqns[i].n+i] = 1.0
            init_sol += [phase1.system.eqns[i].b]

        #Update all n
        for i in range(0,phase1.system.m):
            phase1.system.eqns[i].n += phase1.system.m

        #Set the objective: minimising the new variables, we want them to be 0
        phase1.minmax = "min"
        phase1.objective = [0.0]*(cpy.system.eqns[0].n+phase1.system.m)
        for i in range(0,cpy.system.eqns[0].n):
            for j in range(0,phase1.system.m):
                phase1.objective[i] -= phase1.system.eqns[j].a[i]
            if i < len(cpy.objective):
                cpy.objective[i] -+ phase1.system.eqns[j].a[i]

        #Run the simplex algorithm to find a feasible solution
        phase1wsol = LPwFeasSol(phase1.system, phase1.objective,phase1.minmax,init_sol,cpy.objective)
        phase1wsol.solve()

        #Check we actually found a feasible solution, accounting for some rounding errors
        for i in range(cpy.system.eqns[0].n,phase1wsol.system.eqns[0].n):
            if phase1wsol.sol[i] > 0.001:
                return "No feasible solution"

        #Prepare the new system
        for i in range(0,cpy.system.m):
            cpy.system.eqns[i].a = phase1wsol.system.eqns[i].a[0:cpy.system.eqns[i].n]
            cpy.system.eqns[i].b = phase1wsol.system.eqns[i].b

        #Now run the simplex algorithm seeded with the feasible solution we just found
        pwsol = LPwFeasSol(cpy.system,phase1wsol.secobj,cpy.minmax,phase1wsol.sol)
        r = pwsol.solve()

        if (r == 2):
            return "Unbounded"
        else:
            return [pwsol.sol[i]-pwsol.sol[n+i] for i in range(0,n)]

    #Procedure printmatrix
    #Prints the current matrix, for debugging purposes
    def printmatrix(self):
        print self.objective
        for i in range(0,self.system.m):
            print str(self.system.eqns[i].a) + self.system.eqns[i].op + str(self.system.eqns[i].b)

#Class Linear Program with Feasible Solution, an extension of LinProgram with a solution
#Attributes: those from LinProgram, and:
##sol, a feasible solution (not checked explicitly on instantiation)
##secobj (optional), a secondary objective which is not used for optimisation purposes
##                   but undergoes the same row operations
class LPwFeasSol(LinProgram):

    def __init__(self,system,objective,minmax,sol,secobj=[]):
        LinProgram.__init__(self,system,objective,minmax)
        self.sol = sol
        self.secobj = secobj

    #Procedure onestep
    #Implements one step of the simplex algorithm
    #The program is already assumed to be in normal form
    #We use the Devex pivot rule
    #Output: an integer, either 0 if we are not done iterating, 1 if we are done iterating and
    #        found the optimum, or 2 if we find the problem is unbounded
    def onestep(self):
        n = self.system.eqns[0].n
        m = self.system.m

        #Find the entering variable/pivot column
        min = 0
        entering = -1
        for i in range(0,n):
            if self.sol[i] == 0 and self.objective[i] < 0 and self.objective[i] < min:
                entering = i
                min = self.objective[i]

        #We found the optimum
        if entering == -1:
            return 1

        #Find the leaving variable/pivot row
        leaving = -1
        min = 0
        for i in range(0,m):
            u = self.system.eqns[i].a[entering]
            if u > 0 and (self.system.eqns[i].b/u < min or leaving == -1) and self.system.eqns[i].op == "==":
                leaving = i
                min = self.system.eqns[i].b/u

        #The problem is unbounded
        if leaving == -1:
            return 2

        #Perform the pivot on the leaving row
        d = self.system.eqns[leaving].a[entering]
        self.system.eqns[leaving].a = [x/d for x in self.system.eqns[leaving].a]
        self.system.eqns[leaving].b /= d

        #Perform the pivot on the other rows
        for i in range(0,m):
            if (i != leaving):
                z = self.system.eqns[i].a[entering]
                self.system.eqns[i].a = [self.system.eqns[i].a[j] - z*self.system.eqns[leaving].a[j] for j in range(0,n)]
                self.system.eqns[i].b = self.system.eqns[i].b - z*self.system.eqns[leaving].b

        #Perform the pivot on the objective
        self.objective = [self.objective[j] - self.objective[entering]*self.system.eqns[leaving].a[j] for j in range(0,n)]

        #Perform the pivot on the secondary objective
        self.secobj = [self.secobj[j] - self.secobj[entering]*self.system.eqns[leaving].a[j] for j in range(0,len(self.secobj))]

        #Find the new feasible solution
        self.sol = [0.0]*n

        #Find the new basic variables and set the corresponding values of the solution
        taken = [False]*self.system.m
        for i in range(0,n):
            truth = 0
            k = -1
            for j in range(0,m):
                if self.system.eqns[j].op == "==":
                    if self.system.eqns[j].a[i] != 0.0:
                        truth += 1
                        k = j
                        if self.system.eqns[j].a[i] != 1.0:
                            truth += 1
            if truth == 1 and not taken[k]:
                self.sol[i] = self.system.eqns[k].b
                taken[k] = True

        return 0

    #Procedure solve
    #Performs the simplex algorithm, i.e. iterates self.onestep
    #The program is assumed to be in normal form
    #Output: 2 if the problem is unbounded, and 0 otherwise
    def solve(self):
        r = 0
        while r == 0:
            r = self.onestep()
        return r

    #Procedure evaluate
    #Computes the objective at self.sol
    #Output: the intter product of self.sol and self.objective
    def evaluate(self):
        return LinProgram.evaluate(self,self.sol)

    #Prints the current matrix for debugging purposes
    def printmatrix(self):
        LinProgram.printmatrix(self)
        print self.sol


#The main program code, with an example
ineq1 = LinIneq([3.0,1.0,1.0,4.0,],"<=",12.0)
ineq2 = LinIneq([1.0,-3.0,2.0,3.],"<=",7.0)
ineq3 = LinIneq([2.0,1.0,3.0,-1.],"<=",10.0)
ineq4 = LinIneq([1.0,0.0,0.0,0.0],">=",0.0)
ineq5 = LinIneq([0.0,1.0,0.0,0.0],">=",0.0)
ineq6 = LinIneq([0.0,0.0,1.0,0.0],">=",0.0)
ineq7 = LinIneq([0.0,0.0,0.0,1.0],">=",0.0)

sys1 = LinSystem([ineq1,ineq2,ineq3,ineq4,ineq5,ineq6,ineq7])
program1 = LinProgram(sys1,[2.0,4.0,3.0,1.0],"max")

solution1 = program1.solve()
print "Solution 1: " + str(solution1) + " with value " + str(program1.evaluate(solution1))

#another example

ineq_1 = LinIneq([3.0,2.0,1.0],"==",10.0)
ineq_2 = LinIneq([2.0,5.0,3.0],"==",15.0)
ineq_3 = LinIneq([1.0,0.0,0.0],">=",0.0)
ineq_4 = LinIneq([0.0,1.0,0.0],">=",0.0)
ineq_5 = LinIneq([0.0,0.0,1.0],">=",0.0)

sys2 = LinSystem([ineq_1,ineq_2,ineq_3,ineq_4,ineq_5])
program2 = LinProgram(sys2,[-2.0,-3.0,-4.0],"min")

solution2 = program2.solve()
print "Solution 2: " + str(solution2) + " with value " + str(program2.evaluate(solution2))
