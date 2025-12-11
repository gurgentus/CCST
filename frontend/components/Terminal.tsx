'use client';

import { useEffect, useRef } from 'react';
import { api } from '@/lib/api';

// jQuery and Terminal will be loaded dynamically
declare global {
  interface Window {
    jQuery: any;
    $: any;
  }
}

interface TerminalProps {
  onOutput?: (output: string) => void;
}

export default function Terminal({ onOutput }: TerminalProps) {
  const terminalRef = useRef<HTMLDivElement>(null);
  const termInstanceRef = useRef<any>(null);

  useEffect(() => {
    // Load jQuery and jQuery Terminal dynamically
    const loadScripts = async () => {
      if (typeof window === 'undefined') return;

      // Load jQuery
      if (!window.jQuery) {
        await loadScript('https://code.jquery.com/jquery-3.2.1.min.js');
      }

      // Load jQuery Terminal
      if (!window.jQuery.fn.terminal) {
        await loadScript('https://unpkg.com/jquery.terminal/js/jquery.terminal.min.js');
        await loadStylesheet('https://unpkg.com/jquery.terminal/css/jquery.terminal.min.css');
      }

      // Initialize terminal
      if (terminalRef.current && !termInstanceRef.current) {
        initializeTerminal();
      }
    };

    loadScripts();

    return () => {
      if (termInstanceRef.current) {
        termInstanceRef.current.destroy();
      }
    };
  }, []);

  const loadScript = (src: string): Promise<void> => {
    return new Promise((resolve, reject) => {
      const script = document.createElement('script');
      script.src = src;
      script.onload = () => resolve();
      script.onerror = reject;
      document.head.appendChild(script);
    });
  };

  const loadStylesheet = (href: string): Promise<void> => {
    return new Promise((resolve) => {
      const link = document.createElement('link');
      link.rel = 'stylesheet';
      link.href = href;
      link.onload = () => resolve();
      document.head.appendChild(link);
    });
  };

  const initializeTerminal = () => {
    const $ = window.jQuery;

    const commands = {
      // Matrix operations
      matrix: async function(this: any, name: string, value: string) {
        if (!name || !value) {
          this.error('Usage: matrix <name> <value>');
          return;
        }
        try {
          const res = await api.createMatrix(name, value);
          this.echo(`[[;cyan;]Matrix '${name}' created]`);
          if (onOutput) onOutput(`Matrix '${name}' created`);
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Calculate eigenvectors/eigenvalues
      modes: async function(this: any, mat: string, V: string, E: string) {
        if (!mat || !V || !E) {
          this.error('Usage: modes <matrix> <V_output> <E_output>');
          return;
        }
        try {
          const res = await api.calculateModes(mat, V, E);
          this.echo(`[[;cyan;]Modes calculated: ${V} (eigenvectors), ${E} (eigenvalues)]`);
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Define state-space system
      ss: async function(this: any, A: string, B: string, C: string, D: string, name: string) {
        if (!A || !B || !C || !D || !name) {
          this.error('Usage: ss <A> <B> <C> <D> <name>');
          return;
        }
        try {
          const res = await api.createStateSpace(A, B, C, D, name);
          this.echo(`[[;cyan;]${res.result || 'State-space system created'}]`);
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Define controller
      controller: async function(this: any, A: string, B1: string, B2: string, C: string, D1: string, D2: string, name: string) {
        if (!A || !B1 || !B2 || !C || !D1 || !D2 || !name) {
          this.error('Usage: controller <A> <B1> <B2> <C> <D1> <D2> <name>');
          return;
        }
        try {
          const res = await api.createController(A, B1, B2, C, D1, D2, name);
          this.echo(`[[;cyan;]${res.result || 'Controller created'}]`);
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Feedback system
      feedback: async function(this: any, G: string, K: string, name: string) {
        if (!G || !K || !name) {
          this.error('Usage: feedback <G> <K> <name>');
          return;
        }
        try {
          const res = await api.createFeedback(G, K, name);
          this.echo(`[[;cyan;]${res.result || 'Feedback system created'}]`);
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Designate control input
      control: async function(this: any, index: string, G: string, name: string) {
        if (!index || !G || !name) {
          this.error('Usage: control <index> <G> <name>');
          return;
        }
        try {
          const res = await api.designateControl(parseInt(index), G, name);
          this.echo(`[[;cyan;]${res.result || 'Control input designated'}]`);
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Designate output
      output: async function(this: any, index: string, G: string, name: string) {
        if (!index || !G || !name) {
          this.error('Usage: output <index> <G> <name>');
          return;
        }
        try {
          const res = await api.designateOutput(parseInt(index), G, name);
          this.echo(`[[;cyan;]${res.result || 'Output designated'}]`);
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Display variable graphically
      gdisplay: async function(this: any, name: string) {
        if (!name) {
          this.error('Usage: gdisplay <variable_name>');
          return;
        }
        try {
          const data = await api.getDisplay(name);
          this.echo(`[[;cyan;]Variable '${name}' sent to graphical output below]`);
          if (onOutput) onOutput(JSON.stringify(data));
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Test command
      echo: function(this: any, ...args: string[]) {
        this.echo(args.join(' '));
      },

      // Ping test
      ping: async function(this: any) {
        try {
          const res = await api.ping();
          this.echo(`[[;cyan;]${res.message}]`);
        } catch (e: any) {
          this.error(`Error: ${e.message}`);
        }
      },

      // Help command
      help: function(this: any) {
        const helpText = `
[[;cyan;]Available Commands:]

  [[;yellow;]Matrix Operations:]
    matrix <name> <value>           - Create a matrix
    modes <mat> <V> <E>             - Calculate eigenvalues/eigenvectors
    gdisplay <name>                 - Display variable graphically

  [[;yellow;]Control Systems:]
    ss <A> <B> <C> <D> <name>       - Define state-space system
    controller <A> <B1> <B2> <C> <D1> <D2> <name> - Define controller
    feedback <G> <K> <name>         - Create feedback system
    control <index> <G> <name>      - Designate control input
    output <index> <G> <name>       - Designate output

  [[;yellow;]Utilities:]
    ping                            - Test server connection
    echo <text>                     - Echo text
    help                            - Show this help
    clear                           - Clear terminal

[[;cyan;]Example:]
  > matrix A [[1,2],[3,4]]
  > gdisplay A
`;
        this.echo(helpText);
      }
    };

    termInstanceRef.current = $(terminalRef.current).terminal(commands, {
      greetings: 'Cloud Controls and Simulation Toolbox',
      prompt: '> ',
      height: 300,
      checkArity: false,
    });
  };

  return (
    <div
      ref={terminalRef}
      className="w-full"
    />
  );
}
