import math

""" DATA """

V = [[0.5774, 1.0], [1.0, 1.0], [1.1547, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 0.732], [1.0, 0.1547], [0.732, 0.0], [1.0491, 0.183], [-0.317, 0.549], [1.0, 0.268], [0.183, -0.3169], [0.5491, 1.049], [0.4642, 1.0], [0.0, -0.4226], [0.0, 1.0]]
                                                                                       
EV = [[0, 1], [2, 3], [5, 4], [7, 6], [2, 8], [3, 6], [4, 9], [0, 10], [9, 5], [8, 10], [7, 11], [12, 13], [6, 8], [6, 10], [4, 7], [4, 11], [4, 14], [5, 15], [11, 14], [0, 12], [13, 15], [0, 13], [1, 10], [3, 7], [5, 13]]

""" UTILS """

def sub (v1, v2):
  return [v1[0] - v2[0], v1[1] - v2[1]]


def mod (n, m):
  return ((n % m) + m) % m


""" CYCLES """

def compute_ev_mapping (EV):
  ev_mapping = []
  for i in range(len(EV)):
    ev_mapping.append({
      'color': 0,
      'direction': 0
    })
  return ev_mapping


def compute_angle (P, V):
  point = sub(V, P)
  angle = math.atan2(point[1], point[0])
  return angle


def compute_incidences (V, EV):
  incidences = []
  for i in range(len(V)):
    vertex = V[i]
    incidence = []
    for j in range(len(EV)):
      edge = EV[j]
      endpoint = None
      if (edge[0] == i):
        endpoint = edge[1]
        position = 1
      if (edge[1] == i):
        endpoint = edge[0]
        position = 0
      if (endpoint != None):
          incidence.append({
            'index': j,
            'endpoint': endpoint,
            'angle': compute_angle(vertex, V[endpoint]),
            'edge': edge,
            'position': position
          })
    incidence = sorted(incidence, key=lambda item: item['angle'], reverse=True)
    incidences.append(incidence)
  # print map(lambda l: map(lambda i: i['angle'], l), incidences)
  return incidences

def get_starting_edge (incidences, ev_mapping):
  for e in range(len(ev_mapping)):
    if (ev_mapping[e]['color'] < 2):
      direction = 0 if ev_mapping[e]['direction'] == 1 else 1
      color(ev_mapping, e, direction)
      return {
        'edge': e,
        'position': direction
      }
  return None


def get_next_edge (incidences, edge, position):
  items = incidences[EV[edge][position]]
  for j in range(len(items)):
    item = items[j]
    if (item['index'] == edge):
      out = items[mod(j + 1, len(items))]
      return {
        'edge': out['index'],
        'vertex': out['endpoint'],
        'position': out['position']
      }
  return None

def color (ev_mapping, index, direction):
  ev_mapping[index]['color'] += 1
  ev_mapping[index]['direction'] = 1 if direction == 1 else -1


def find_cycles (V, EV):
  ev_mapping = compute_ev_mapping(EV)
  incidences = compute_incidences(V, EV)
  V_cycles = []
  E_cycles = []
  counter = 1
  start = get_starting_edge(incidences, ev_mapping)
  while (start != None):
    V_cycle = [EV[start['edge']][mod(start['position'] + 1, 2)], EV[start['edge']][start['position']]]
    E_cycle = [start['edge']]
    next = get_next_edge(incidences, start['edge'], start['position'])
    while (next['edge'] != start['edge']):
      V_cycle.append(next['vertex'])
      E_cycle.append(next['edge'])
      color(ev_mapping, next['edge'], next['position'])
      next = get_next_edge(incidences, next['edge'], next['position'])
    E_cycles.append(E_cycle)
    V_cycles.append(V_cycle)
    print '############## CYCLE ', counter
    print 'EDGES:', E_cycle
    print 'VERTICES:', V_cycle
    print 'START', 'edge:', start['edge'], 'position:', start['position']
    print 'COUNTER:', map(lambda e: e['color'], ev_mapping)
    print '\n'
    start = get_starting_edge(incidences, ev_mapping)
    counter += 1
  return {
    'v_cycles': V_cycles,
    'e_cycles': E_cycles,
    'ev_mapping': ev_mapping
  }


""" MAIN """

cycles_data = find_cycles(V, EV)
print '############## OUTPUT'
print 'EDGES:'
print cycles_data['e_cycles']
print '\n'
print 'VERTICES:'
print cycles_data['v_cycles']


