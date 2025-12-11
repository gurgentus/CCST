// API utility for communicating with FastAPI backend

const API_BASE = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

export interface ApiResponse<T = any> {
  result?: string;
  message?: string;
  detail?: string;
  [key: string]: any;
}

export class ApiClient {
  private baseUrl: string;

  constructor(baseUrl: string = API_BASE) {
    this.baseUrl = baseUrl;
  }

  private async request<T>(
    endpoint: string,
    options: RequestInit = {}
  ): Promise<T> {
    const url = `${this.baseUrl}${endpoint}`;

    try {
      const response = await fetch(url, {
        ...options,
        headers: {
          'Content-Type': 'application/json',
          ...options.headers,
        },
      });

      if (!response.ok) {
        const error = await response.json().catch(() => ({
          detail: 'Request failed'
        }));
        throw new Error(error.detail || `HTTP ${response.status}`);
      }

      return await response.json();
    } catch (error) {
      console.error('API request failed:', error);
      throw error;
    }
  }

  // Control API endpoints
  async createMatrix(name: string, value: string): Promise<ApiResponse> {
    return this.request('/api/control/matrix', {
      method: 'POST',
      body: JSON.stringify({ name, value }),
    });
  }

  async calculateModes(mat: string, V: string, E: string): Promise<ApiResponse> {
    return this.request('/api/control/modes', {
      method: 'POST',
      body: JSON.stringify({ mat, V, E }),
    });
  }

  async createStateSpace(
    A: string,
    B: string,
    C: string,
    D: string,
    name: string
  ): Promise<ApiResponse> {
    return this.request('/api/control/ss', {
      method: 'POST',
      body: JSON.stringify({ A, B, C, D, name }),
    });
  }

  async createController(
    A: string,
    B1: string,
    B2: string,
    C: string,
    D1: string,
    D2: string,
    name: string
  ): Promise<ApiResponse> {
    return this.request('/api/control/controller', {
      method: 'POST',
      body: JSON.stringify({ A, B1, B2, C, D1, D2, name }),
    });
  }

  async createFeedback(G: string, K: string, name: string): Promise<ApiResponse> {
    return this.request('/api/control/feedback', {
      method: 'POST',
      body: JSON.stringify({ G, K, name }),
    });
  }

  async designateControl(index: number, G: string, name: string): Promise<ApiResponse> {
    return this.request('/api/control/control', {
      method: 'POST',
      body: JSON.stringify({ index, G, name }),
    });
  }

  async designateOutput(index: number, G: string, name: string): Promise<ApiResponse> {
    return this.request('/api/control/output', {
      method: 'POST',
      body: JSON.stringify({ index, G, name }),
    });
  }

  // Display API
  async getDisplay(name: string): Promise<any> {
    return this.request(`/gdisplay?name=${encodeURIComponent(name)}`);
  }

  // Test endpoint
  async ping(): Promise<ApiResponse> {
    return this.request('/api/ping');
  }
}

export const api = new ApiClient();
