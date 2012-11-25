import strife
import scipy as sp
import sys
import os
import signal
import time

# TODO: Handler of signals.
def handler_maker(a_game):
  def handler(signum, frame):
      print 'Signal handler called with signal', signum
      a_game.save_h5()
      print 'game saved'
      raise
  return handler

def go(a):
    signal.signal(signal.SIGINT, handler_maker(a))
    t = time.time()
    every = 30*60
    # TODO: Maybe add f and d somehow like in printf? {0}f {1}d
    print "t: {0:f}, steps thus far: {1}".format(t, a.step_count)
    steps_a = a.step_count
    while a.step_count <= a.generations:
        a.nextstep()
        delta_t = time.time() - t
        if delta_t > every:
            t = time.time()
            a.save_h5()
            steps_delta = a.step_count - steps_a
            steps_a = a.step_count
            eta = 1.0 * delta_t / (steps_delta+1) * (a.steps_final - a.step_count)
            print "t: {0:f}, approx. time to fin: {1:f}".format(t, eta)
            print "steps taken = {0}, steps since last save = {1}".format(a.step_count, steps_delta)
            sys.exit(1)

if __name__ == '__main__':
      config = strife.default_config
      for i, arg in enumerate(sys.argv):
          if arg in ('--datafile', '-d'):
              config['data_filename'] = sys.argv[i+1]
          if arg in ('--config', '-c'):
              conf_filename = sys.argv[i+1]
              config = strife.load_config(conf_filename)
      if os.path.exists(config['data_filename']):
          game = strife.gameOfStrife()
          game.load_h5()
          go(game)
      else:
          game = strife.gameOfStrife(config)
          game.save_h5()
          go(game)
      game.save_h5()
      sys.exit(0)
