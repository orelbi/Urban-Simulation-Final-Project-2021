"""תרגיל מסכם בקורס סימולציה עירונית"""

# Imports relevant libraries
from random import choice, random, randint
import shapefile
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# class for creating shop objects
class Shop:
    # לכל חנות יש נתונים: מיקום, עלות תפעול, איכות מוצר, מחיר
    # attributes
    location = None  # tuple - x,y
    operation_costs = None  # costs for running the shop
    quality = None  # quality of the products sold in the shop
    price = None  # price of products sold in the shop

    def __init__(self, location, quality, price, operation_costs):
        # לכל חנות יש את הנתונים שצויינו קודם ונתון בדבר הרווחים
        self.location = location
        self.quality = quality
        self.price = price
        self.operation_costs = operation_costs
        self.revenue = 0

    def dist(self, l):
        # compute aerial distance between object and location l
        # חישוב מרחק אווירי בין החנות לבין מיקום הסוכן
        return ((self.location[0] - l[0]) ** 2 + (self.location[1] - l[1]) ** 2) ** 0.5

    def step(self):
        # חישוב רווח כללי בסימולציה
        model.revenues += self.revenue  # increase model variable to document revenue

        # אם הרווח קטן מעלות התפעול החנות תוריד מחיר
        if self.revenue < self.operation_costs:  # if the shop is not profitable, decrease price
            self.price -= 1  # decrease price
            # אם המחיר 0 משמע החנות לא תהיה יעילה ולכן עליה להיסגר
            if self.price == 0:  # if price is not positive - remove from model and agents' maps
                model.venues.remove(self)
                for a in model.agents:
                    if self in a.cognitive_map:
                        del a.cognitive_map[self]
        # אם הרווחים גדולים פי 2 ומעלה מעלויות התפעול החנות תעלה את המחיר
        # כדי למקסם רווחים
        elif self.revenue >= 2 * self.operation_costs:  # if profits are high - increase price
            self.price += 1
        # איפוס רווחים לקראת הצעד הבא במודל
        self.revenue = 0  # ireset revenue for next round


# class for creating individual agent objects
class Agent:
    # לכל סוכן יש נתונים של מיקום, תקציב לנסיעה ולקנייה, כמות הצריכה הרצויה
    # attributes
    location = None  # location of agent - tuple of (x,y)
    budget = None  # budget for travel and purchasing products
    consumption = None  # desired amount of product to be purchased

    # auxilliary variable that will contain the agent's knowledge about the world
    # משתנה עזר שיכיל את הידע של הסוכן על העולם
    cognitive_map = None

    def __init__(self, location, budget, consumption):
        self.location = location
        self.budget = budget
        self.consumption = consumption

    def compute_utility(self, b, quality):  # input - a shop, quality level
        # מחשב את התועלת, הקלט זה חנות B ואיכות מוצר
        # compute cost of travelling to shop
        # חישוב עלות נסיעה לחנות
        transport_costs = b.dist(self.location) * model.transportation_costs
        # compute the amount the agent will purchase given transport costs, budget, and price
        # חישוב מה הכמות המקסימלית של המוצר שהסוכן יכול לרכוש
        consumption_amount = int((self.budget - transport_costs) / b.price)
        if consumption_amount > self.consumption:  # make sure consumption amount does not exceed the desired amount
            consumption_amount = self.consumption
        if consumption_amount < 0:  # make sure consumption amount is not negative
            # אם הערך של הכמות שלילית הסוכן לא ירכוש כלל את המוצר
            consumption_amount = 0
        # חישוב תועלת
        # TODO 1
        utility = ((consumption_amount/self.consumption)*(1-model.w_quality)) + quality*model.w_quality
        # החזרת התועלת
        return utility

    def step(self):
        # find known shops that are believed to produce enough utility
        # חישוב התועלת מבוסס על חנות ועל איכות המוצר בעיניו של הסוכן
        # קודם עוברים על חנויות המצויות במפה הקוגנטיבית
        candidates = [b for b in self.cognitive_map if
                      self.compute_utility(b, self.cognitive_map[b]) >= model.utility_threshold]
        if len(candidates) > 0:  # if found such shops - randomly draw one
            b = choice(candidates)
        else:
            # אם החיפוש העלה חרס נעבור על חנויות פתוחות שאינן במפה הקוגנטיבית
            # find shops not known to the agent
            candidates = [b for b in model.venues if b not in self.cognitive_map]
            if len(candidates) > 0:  # if there are such shops - draw one randomly
                b = choice(candidates)
                # נכניס את החנות למפה הקוגנטיבית וניתן ערך רנדומלי לאיכות
                # לפי דעתו של הסוכן
                self.cognitive_map[b] = random()  # add to cognitive map and randomly set beliefs regarding quality
            else:  # otherwise - exit function
                return
        # compute transport costs and the amount purchased by the agent at the shop
        # חישוב עלות נסיעה
        transport_costs = b.dist(self.location) * model.transportation_costs
        # חישוב כמות המוצר שסוכן יכול לקנות
        consumption_amount = int((self.budget - transport_costs) / b.price)
        # סידור כמות המוצר שאדם קונה
        # make sure amount is positive and not bigger than desired amount
        if consumption_amount < 0:
            consumption_amount = 0
        if consumption_amount > self.consumption:
            consumption_amount = self.consumption
        # עדכון הרווח של החנות באותו צעד
        b.revenue += consumption_amount * b.price  # update shop's revenue with the consumption amount * price
        # TODO 2
        # עדכון איכות המוצר להיות ממוצע האמונה באיכות
        # ומהאיכות המוגדרת במודל למוצר
        new_belief = (b.quality + self.cognitive_map[b])/2
        self.cognitive_map[b] = new_belief
        # עדכון סכום תועלות על בסיס הוספת תועלת מבוססת חנות ואיכות מוצר ע"פ
        #  המודל
        model.utilities += self.compute_utility(b, b.quality)  # get actual utility from the purchase
        # הגדלת כמות הביקורים ב-1
        model.visits += 1


# class for creating model objects
class Model:
    # הגדרת מאפיינים: כמות צריכה רצויה, משקל איכות המוצר, תועלת מינימלית
    # שסוכן ירצה להפיק, רשימת חנויות פתוחות, סוכנים, תועלות, רווחים,
    # משתנה עזר לשמירה על תוצאות, מספר סימולציה, מספר ביקורים בחנויות
    transportation_costs = None  # model parameter - cost per distance unit
    w_quality = None  # weight given by agents to quality of products
    utility_threshold = None  # minimal utility agents are requiring from shops
    venues = None  # list of shops
    agents = None  # list of agents
    utilities = None  # auxilliary variable used to document the sum of agent's utilities
    revenues = None  # auxilliary variable used to document the sum of shop revenues
    data = None  # auxilliary variable used to save outputs
    sim = None  # documents the number of the simulation
    visits = None  # auxilliary variable used to document the number of time agents visited shops

    def __init__(self, transportation_costs, shape, sim):
        self.transportation_costs = transportation_costs
        self.w_quality = 0.5
        self.utility_threshold = 0.6
        self.sim = sim
        self.venues = []
        self.agents = []
        self.data = []

        shp = shapefile.Reader(shape)  # read shapefile
        for i in range(len(shp.shapes())):
            if shp.record(i)['landuse'] == 'commercial':  # if record is commercial, create a shop
                self.venues.append(Shop(shp.shape(i).points[0],  # location
                                        random(),  # quality
                                        randint(10, 50),  # price
                                        randint(1000, 10000)))  # operation costs
            else:  # if record is residential, create an agent
                self.agents.append(Agent(shp.shape(i).points[0],  # location
                                         randint(200, 400),  # budget
                                         randint(1, 10)))  # desired consumption amount

        for a in self.agents:
            # לכל סוכן יש במפה הקוגנטיבית 3 חנויות
            # randomly add 3 shops to the cognitive map and set a random beliefs regarding quality
            a.cognitive_map = {choice(self.venues): random() for i in range(3)}

    def run_model(self):
        for i in range(100):
            self.utilities = 0
            self.revenues = 0
            self.visits = 0
            for a in self.agents:  # agents consume
                a.step()
            for b in self.venues:  # shops update prices
                b.step()

            # append to data - transportation costs, simulation number, step number,
            # number of open shops, number of agent visits, average number of shops known to agents,
            # average utility derived from purchase per agent, average revenue per shop, average price per shop
            # נוסיף את הנתונים הבאים: עלויות נסיעה, מספר סימולציה,
            # מספר צעדים, מספר חנויות פתוחות, מספר ביקורים, מספר חנויות
            # ממוצע ידוע לסוכנים, תועלת ממוצעת לכל סוכן מרכישה, רווח ממוצע
            # לכל חנות, מחיר ממוצע לחנות

            self.data.append([self.transportation_costs, self.sim, i, len(self.venues), self.visits,
                              round(sum([len(a.cognitive_map) for a in self.agents]) / len(self.agents), 2),
                              round(self.utilities / len(self.agents), 2),
                              round(self.revenues / len(self.venues), 2),
                              round(sum([b.price for b in self.venues]) / len(self.venues), 2)])


results = []
for j in range(0, 10, 2):  # iterate over values of transportation costs
    for i in range(10):  # simulate 10 times
        print(j / 10, i)
        model = Model(j / 10, 'bldgs_points', i)
        model.run_model()
        for d in model.data:  # add outputs to results
            results.append(d)
# transform results to data frame
df = pd.DataFrame(results, columns=['Transport costs', 'Simulation', 'Step', 'Venues',
                                    'Visits', 'Known venues', 'Mean utility', 'Mean revenue',
                                    'Mean price'])
# TODO 3
# יצירת פלט של line plot עבור כל עמודה בטבלה
# create line plots for each of the documented variables

for c in ['Venues', 'Visits', 'Known venues', 'Mean utility', 'Mean revenue',
          'Mean price']:

    fig, axs = plt.subplots(1, 1)

    sns.lineplot(data=df, x="Step", y=c, hue='Transport costs',
                 legend='full', palette='tab10')

    axs.set_title("Line Plot " + str(c))
    fig.tight_layout()
    fig.savefig(str(c) + ".png")
