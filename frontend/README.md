# Kinection — Frontend

Web frontend for Kinection, a personal ancient-DNA comparison platform.
Vite + React + TypeScript SPA, styled with Tailwind and shadcn/ui-style
components. Targets Cloudflare Pages.

## Run

```bash
cd frontend
npm install
npm run dev        # http://localhost:5173
```

Build for production:

```bash
npm run build      # tsc -b && vite build -> dist/
npm run preview    # serve the production build locally
```

## Environment

The API base URL is read from `VITE_API_BASE_URL`. Copy `.env.example` to
`.env` and adjust as needed:

```
VITE_API_BASE_URL=http://localhost:8787
```

If unset, the client defaults to `http://localhost:8787` (local Worker dev).

## Flow

1. `/` — pick a raw DNA `.txt` file, optional label, upload. The client calls
   `POST /uploads/url` for a presigned R2 URL, `PUT`s the file directly to R2,
   then calls `markUploadComplete()`.
2. `/jobs/:id` — polls `GET /jobs/:id` every 3s; shows status; has a
   "Delete my upload" button (`DELETE /jobs/:id/upload`).
3. `/jobs/:id/report` — fetches `report.json` + `map_data.geojson` and renders
   the report sections.

## Known backend follow-ups

- **`POST /jobs/:id/upload-complete` does not exist yet.** After the file `PUT`
  succeeds, the job stays in status `uploading`. `api.markUploadComplete()`
  calls this route and swallows a 404 as a known-pending gap (logged as a
  warning). Until the backend implements it, the status page will sit on
  `uploading`. See the `TODO(backend)` comment in `src/api/client.ts`.

## Placeholders (intentionally incomplete)

- **Map** (`MapPlaceholder.tsx`) — renders a react-leaflet map with OSM tiles
  and loads `map_data.geojson`, but does not yet plot features as markers.
- **PCA** (`PcaPlaceholder.tsx`) — static "coming soon" card; PCA output is
  currently empty in reports.

## Deploy (Cloudflare Pages)

Build command `npm run build`, output directory `dist`. Set
`VITE_API_BASE_URL` to the deployed Worker URL as a Pages build-time env var.
This is an SPA, so configure a catch-all rewrite to `/index.html` for client
routing.
