// Terminal command types

export interface TerminalCommand {
  name: string;
  description: string;
  usage: string;
  handler: (term: any, ...args: string[]) => Promise<void>;
}

export interface DisplayData {
  what: 'matrix' | 'color_matrix' | 'ss' | 'orbit' | 'controller' | 'feedback' | 'plot' | 'image';
  value?: any;
  A?: { value: number[][] };
  B?: { value: number[][] };
  C?: { value: number[][] };
  D?: { value: number[][] };
  B1?: { value: number[][] };
  B2?: { value: number[][] };
  D1?: { value: number[][] };
  D2?: { value: number[][] };
  h?: number;
  a?: number;
  e?: number;
  Omega?: number;
  i?: number;
  omega?: number;
  head?: string;
  script?: string;
  div?: string;
}

export interface TerminalOutput {
  type: 'success' | 'error' | 'info';
  message: string;
}
